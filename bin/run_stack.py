#!/usr/bin/env python

import sys, os
import datetime
import shutil
import optparse
import collections
import multiprocessing
import hsc.meas.mosaic.mosaicLib  as hscMosaic
import hsc.meas.mosaic.stack             as stack

import lsst.obs.hscSim as obsHsc
import lsst.obs.suprimecam              as obsSc
import lsst.pipette.readwrite as pipReadWrite
import lsst.pipette.runHsc as runHsc
import lsst.afw.image as afwImage
try:
    from IPython.core.debugger import Tracer
    debug_here = Tracer()
except:
    pass


WarpInputs = collections.namedtuple('WarpInputs', ['fileIO', 'f', 'skipMosaic', 'instrument', 'rerun'])

def runStackWarp(warpInputs):

    fileIO     = warpInputs.fileIO
    f          = warpInputs.f
    skipMosaic = warpInputs.skipMosaic
    instrument = warpInputs.instrument
    rerun      = warpInputs.rerun
    
    wcsDic, dims, fscale = stack.readParamsFromFileList([f], skipMosaic=skipMosaic)

    if instrument.lower() in ["hsc"]:
        mapper = obsHsc.HscSimMapper(rerun=rerun)
    elif instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
        mapper = obsSc.SuprimecamMapper(rerun=rerun)
    ioMgr = pipReadWrite.ReadWrite(mapper, ['visit', 'ccd'], config={})
    
    trueSigma = -1.0
    try:
	warpResult = stack.stackMeasureWarpedPsf(f, wcsDic[0], ioMgr=ioMgr,
					    skipMosaic=skipMosaic, fileIO=fileIO)
	if fileIO:
	    trueSigma = warpResult
	else:
	    psf, trueSigma = warpResult
    finally:
	return trueSigma


Inputs = collections.namedtuple('Inputs', ['rerun', 'instrument', 'ix', 'iy', 'subImgSize', 'stackId', 'imgMargin', 'fileIO', 'workDir', 'skipMosaic', 'filter', 'matchPsf'])

def runStackExec(inputs):
    rerun = inputs.rerun
    instrument = inputs.instrument
    ix = inputs.ix
    iy = inputs.iy
    subImgSize = inputs.subImgSize
    stackId = inputs.stackId
    imgMargin = inputs.imgMargin
    fileIO = inputs.fileIO
    workDir = inputs.workDir
    skipMosaic = inputs.skipMosaic
    filter = inputs.filter
    matchPsf= inputs.matchPsf

    print 'runStackExec ', ix, iy
    if instrument.lower() in ["hsc"]:
        mapper = obsHsc.HscSimMapper(rerun=rerun)
    elif instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
        mapper = obsSc.SuprimecamMapper(rerun=rerun)
    ioMgr = pipReadWrite.ReadWrite(mapper, ['visit', 'ccd'], config={})

    try:	
        stack.stackExec(ioMgr, ix, iy, stackId,
                        subImgSize, imgMargin,
                        fileIO=fileIO,
                        workDir=workDir,
                        skipMosaic=skipMosaic,
                        filter=filter, matchPsf=matchPsf)
    finally:
        return

def main():
    parser = optparse.OptionParser()
    parser.add_option("-r", "--rerun",
                      type=str, default=None,
                      help="rerun name to take corrected frames from and write stack images to.")
    parser.add_option("-I", "--instrument",
                      type=str, default='suprimecam',
                      help="instument to treat (hsc or suprimecam)")
    parser.add_option("-p", "--program",
                      type=str, default=None,
                      help="program name (e.g. COSMOS_0)")
    parser.add_option("-f", "--filter",
                      type=str, default=None,
                      help="filter name (e.g. W-S-I+)")
    parser.add_option("-d", "--dateObs",
                      type=str, default=None,
                      help="(optional) dataObs (e.g. 2008-11-27)")
    parser.add_option("-t", "--threads",
                      type=int, default=1,
                      help="(optional) Number of threads")
    parser.add_option("-w", "--workDir",
                      type=str, default=None,
                      help="working directory to store files (e.g. /data/yasuda/DC2/sc)")
    parser.add_option("-W", "--workDirRoot",
                      type=str, default='.',
                      help="(optional) root working directory (working dir will be root/program/filter")
    parser.add_option("-s", "--destWcs",
                      type=str, default=None,
                      help="destination wcs")
    parser.add_option("--pScale",
                      type="float", default=0.0,
                      help="destination pixel scale in arcsec")
    parser.add_option("-m", "--doMatchPsf",
		      default=False, action='store_true',
		      help="match PSFs before stacking (default=%default)")
    
    (opts, args) = parser.parse_args()

    if not opts.rerun or not opts.program or not opts.filter:
        parser.print_help()
        raise SystemExit("failed to parse arguments")

    sys.argv = [sys.argv[0]] + args
    print "rerun=%s, instrument=%s, program=%s, filter=%s, dateObs=%s, destWcs=%s, pScale=%f, threads=%d, args=%s " % \
        (opts.rerun, opts.instrument, opts.program, opts.filter, opts.dateObs, opts.destWcs, opts.pScale, opts.threads, sys.argv)

    run(rerun=opts.rerun, instrument=opts.instrument, program=opts.program,
        filter=opts.filter, dateObs=opts.dateObs, destWcs=opts.destWcs,
        pScale=opts.pScale,
        workDir=opts.workDir, workDirRoot=opts.workDirRoot, threads=opts.threads, doMatchPsf=opts.doMatchPsf)
    
def run(rerun=None, instrument=None, program=None, filter=None, dateObs=None, 
        destWcs=None, pScale=0.0, workDir=None, workDirRoot=None, threads=None, doMatchPsf=False):
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    if instrument.lower() in ["hsc"]:
        mapper = obsHsc.HscSimMapper(rerun=rerun)
        ccdIds = range(100)
    elif instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
        mapper = obsSc.SuprimecamMapper(rerun=rerun)
        ccdIds = range(10)

    config = {}
    ioMgr = pipReadWrite.ReadWrite(mapper, ['visit', 'ccd'], config=config)
    if (dateObs == None):
        frameIds = ioMgr.inButler.queryMetadata('calexp', None, 'visit', dict(field=program, filter=filter))
        pointings = ioMgr.inButler.queryMetadata('calexp', None, 'pointing', dict(field=program, filter=filter))
    else:
        frameIds = ioMgr.inButler.queryMetadata('calexp', None, 'visit', dict(field=program, filter=filter, dateObs=dateObs))
        pointings = ioMgr.inButler.queryMetadata('calexp', None, 'pointing', dict(field=program, filter=filter, dateObs=dateObs))
    print frameIds
    print pointings

    if workDirRoot:
        workDir = os.path.join(workDirRoot, program, filter)
    subImgSize = 4096
    imgMargin = 256
    fileIO = True
    writePBSScript = True
    skipMosaic = False
    stackId = pointings[0]

    if (len(sys.argv) == 1):
        fileList = []
        for frameId in frameIds:
            for ccdId in ccdIds:
                try:
                    fname = ioMgr.read('calexp_filename', dict(visit=frameId, ccd=ccdId))[0][0]
                except Exception, e:
                    print "failed to get file for %s:%s" % (frameId, ccdId)
                    continue
                if os.path.isfile(fname):
                    fileList.append(fname)
                else:
                    print "file %s does not exist " % (fname)

#        stack.stack(ioMgr, fileList, subImgSize, stackId, imgMargin=256, fileIO=True,
#                    workDir=workDir, skipMosaic=skipMosaic, filter=filter,
#                    destWcs=destWcs)
        try:
            os.makedirs(workDir)
        except OSError:
            print "Working directory already exists"

        if destWcs != None:
            destWcs = os.path.abspath(destWcs)

	#############################################################
	#  Init
	#############################################################
        nx, ny, fileList, wcs = stack.stackInit(ioMgr, fileList, subImgSize, imgMargin,
                                 fileIO,
                                 workDir=workDir,
                                 skipMosaic=skipMosaic,
                                 destWcs=destWcs)


	#############################################################
	#  measure warped PSF
	#############################################################
	ixs = range(nx)
	iys = range(ny)

	maxWidth = (max(ixs)+1)*subImgSize
	maxHeight = (max(iys)+1)*subImgSize

	
	wcsDic, dims, fscale = stack.readParamsFromFileList(fileList,
						      skipMosaic=skipMosaic)

	# get a list of files which overlap the ixs,iys requested
	# This makes it possible to hack the code easily to run a single ix,iy for debugging
	fileList, wcsDic, dims, fscale = \
		  stack.cullFileList(fileList, wcsDic, ixs, iys, wcs, subImgSize, maxWidth, maxHeight, dims, fscale, nx, ny)

	warpInputs = list()
	for f in fileList:
	    warpInputs.append(WarpInputs(fileIO=fileIO, f=f,
					 skipMosaic=skipMosaic, instrument=instrument, rerun=rerun))

	# process the job
        pool = multiprocessing.Pool(processes=threads)
        sigmas = pool.map(runStackWarp, warpInputs)
        pool.close()
        pool.join()

	# use the largest sigma to determine the size of the double Gaussian PSF we want to match to.
	# we must *degrade* to worst seeing
	print "sigmas: ", sigmas
	matchPsf = None
        if sigmas:
            maxSigma = max(sigmas)
            sigma1 = maxSigma
            sigma2 = 2.0*maxSigma
            kwid = int(4.0*sigma2) + 1
            peakRatio = 0.1
            matchPsf = ['DoubleGaussian', kwid, kwid, sigma1, sigma2, peakRatio]


	
	#############################################################
	#  Exec
	#############################################################
        inputs = list()

	for iy in iys:
	    for ix in ixs: 
                inputs.append(Inputs(rerun=rerun,instrument=instrument,ix=ix, iy=iy,
                                     stackId=stackId, subImgSize=subImgSize, imgMargin=imgMargin,
                                     fileIO=fileIO, workDir=workDir, skipMosaic=skipMosaic, filter=filter,
				     matchPsf=matchPsf
				     ))

        pool = multiprocessing.Pool(processes=threads)
        pool.map(runStackExec, inputs)
        pool.close()
        pool.join()

	
	#############################################################
	#  Stack *END*
	#############################################################
        expStack = stack.stackEnd(ioMgr, stackId, subImgSize, imgMargin,
                                  fileIO=fileIO, width=maxWidth, height=maxHeight,
				  nx=max(ixs)+1, ny=max(iys)+1,
                                  workDir=workDir, filter=filter)

	expStack.writeFits('expStack.fits')
	
    elif (len(sys.argv) == 3):
        ix = int(sys.argv[1])
        iy = int(sys.argv[2])

        stack.stackExec(ioMgr, ix, iy, stackId, subImgSize, imgMargin,
                        fileIO=fileIO,
                        workDir=workDir, skipMosaic=skipMosaic,
                        filter=filter)

    else:
        if (sys.argv[1] == "Init"):
            fileList = []
            for frameId in frameIds:
                for ccdId in ccdIds:
                    try:
                        fname = ioMgr.read('calexp_filename', dict(visit=frameId, ccd=ccdId))[0][0]
                    except Exception, e:
                        print "failed to get file for %s:%s" % (frameId, ccdId)
                        continue
                    #fname = mgr.getCorrFilename(int(frameId), int(ccdId))
                    if os.path.isfile(fname):
                        fileList.append(fname)
                    else:
                        print "file %s does not exist " % (fname)
                    
            try:
                os.makedirs(workDir)
            except OSError:
                print "Working directory already exists"
            #productDir = os.environ.get("hscMosaic".upper() + "_DIR", None)
            #shutil.copyfile(os.path.join(productDir, "example/run_stack.py"),
            #                os.path.join(workDir, "run_stack.py"))
            
            stack.stackInit(ioMgr, fileList, subImgSize, imgMargin, 
                            fileIO, writePBSScript,
                            workDir=workDir, skipMosaic=skipMosaic,
                            rerun=rerun, instrument=instrument, program=program,
                            filter=filter, dateObs=dateObs, destWcs=destWcs,
                            pScale=pScale)

        elif (sys.argv[1] == "End"):
            stack.stackEnd(ioMgr, stackId, subImgSize, imgMargin,
                           fileIO=fileIO,
                           workDir=workDir, filter=filter)
	    
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

if __name__ == '__main__':
    main()
    