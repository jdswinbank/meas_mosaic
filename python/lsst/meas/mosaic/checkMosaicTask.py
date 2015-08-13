import os
import math
import numpy

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import lsst.pex.config                  as pexConfig
import lsst.pipe.base                   as pipeBase
from lsst.meas.mosaic.mosaicTask import MosaicTask
from lsst.meas.mosaic.mosaicTask import MosaicConfig
from lsst.meas.mosaic.mosaicTask import MosaicRunner
import lsst.meas.mosaic.mosaicLib       as measMosaic
import lsst.meas.astrom                 as measAstrom
import lsst.afw.image                   as afwImage
import lsst.afw.geom                    as afwGeom
import lsst.afw.table                   as afwTable

from lsst.meas.photocal.colorterms import ColortermLibraryConfig

class CheckMosaicConfig(MosaicConfig):
    maxMag = pexConfig.Field(
        doc="Maximum magnitude for delta mag histogram",
        dtype=float,
        default=21.0)

class CheckMosaicTask(MosaicTask):

    RunnerClass = MosaicRunner
    canMultiprocess = False
    ConfigClass = CheckMosaicConfig
    _DefaultName = "CheckMosaic"

    def _getConfigName(self):
        return None

    def makeDiffPosFlux(self, allMat, allSource, wcsDic, calibDic, ffpDic):
        dx_m = list()
        dy_m = list()
        dx_s = list()
        dy_s = list()
        m0_m = list()
        dm_m  = list()
        m0_s = list()
        dm_s  = list()
        for ss in allMat:
            ra_cat = ss[0].getRa().asDegrees()
            dec_cat = ss[0].getDec().asDegrees()
            mag_cat = -2.5*math.log10(ss[0].getFlux())
            for j in range(1,len(ss)):
                if (ss[j].getFlux() > 0 and ss[j].getFlux() != float('inf')):
                    iexp = ss[j].getExp()
                    ichip = ss[j].getChip()
                    x, y = ss[j].getX(), ss[j].getY()
                    ra, dec = wcsDic[iexp][ichip].pixelToSky(x, y).getPosition()
                    dx_m.append((ra - ra_cat) * 3600)
                    dy_m.append((dec - dec_cat) * 3600)
                    mag = 2.5*math.log10(calibDic[iexp][ichip].getFluxMag0()[0]/ss[j].getFlux())
                    mcor = ffpDic[iexp][ichip].eval(x,y)
                    jcor = -2.5*math.log10(wcsDic[iexp][ichip].pixArea(afwGeom.Point2D(x, y)) / wcsDic[iexp][ichip].pixelScale().asDegrees()**2)
                    m0_m.append(mag_cat)
                    dm_m.append(mag+mcor+jcor-mag_cat)

            if len(ss) > 2:
                n = 0
                ra_cat = 0.0
                dec_cat = 0.0
                S  = 0.0
                Sx = 0.0
                ra_source = list()
                dec_source = list()
                mag_source = list()
                for j in range(1,len(ss)):
                    if (ss[j].getFlux() > 0 and ss[j].getFlux() != float('inf')):
                        iexp = ss[j].getExp()
                        ichip = ss[j].getChip()
                        x, y = ss[j].getX(), ss[j].getY()
                        ra, dec = wcsDic[iexp][ichip].pixelToSky(x, y).getPosition()
                        ra_source.append(ra)
                        dec_source.append(dec)
                        n += 1
                        ra_cat += ra
                        dec_cat += dec

                        mag = 2.5*math.log10(calibDic[iexp][ichip].getFluxMag0()[0]/ss[j].getFlux())
                        err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                        mcor = ffpDic[iexp][ichip].eval(x,y)
                        jcor = -2.5*math.log10(wcsDic[iexp][ichip].pixArea(afwGeom.Point2D(x, y)) / wcsDic[iexp][ichip].pixelScale().asDegrees()**2)
                        mag_source.append(mag+mcor+jcor)
                        Sx += (mag+mcor+jcor) / (err*err)
                        S  += 1. / (err*err)

                if n != 0 and S != 0:
                    ra_cat /= n
                    dec_cat /= n
                    mag_cat = Sx / S
                    for ra, dec, mag in zip(ra_source, dec_source, mag_source):
                        dx_s.append((ra - ra_cat) * 3600)
                        dy_s.append((dec - dec_cat) * 3600)
                        m0_s.append(mag_cat)
                        dm_s.append(mag - mag_cat)

        dx_m = numpy.array(dx_m)
        dy_m = numpy.array(dy_m)
        m0_m = numpy.array(m0_m)
        dm_m = numpy.array(dm_m)

        for ss in allSource:
            n = 0
            ra_cat = 0.0
            dec_cat = 0.0
            Sx = 0.0
            S  = 0.0
            ra_source = list()
            dec_source = list()
            mag_source = list()
            for j in range(1,len(ss)):
                if (ss[j].getFlux() > 0 and ss[j].getFlux() != float('inf')):
                    iexp = ss[j].getExp()
                    ichip = ss[j].getChip()
                    x, y = ss[j].getX(), ss[j].getY()
                    ra, dec = wcsDic[iexp][ichip].pixelToSky(x, y).getPosition()
                    ra_source.append(ra)
                    dec_source.append(dec)
                    n += 1
                    ra_cat += ra
                    dec_cat += dec

                    mag = 2.5*math.log10(calibDic[iexp][ichip].getFluxMag0()[0]/ss[j].getFlux())
                    err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                    mcor = ffpDic[iexp][ichip].eval(x,y)
                    jcor = -2.5*math.log10(wcsDic[iexp][ichip].pixArea(afwGeom.Point2D(x, y)) / wcsDic[iexp][ichip].pixelScale().asDegrees()**2)
                    mag_source.append(mag+mcor+jcor)
                    Sx += (mag+mcor+jcor) / (err*err)
                    S  += 1. / (err*err)

            if n != 0:
                ra_cat /= n
                dec_cat /= n
                mag_cat = Sx / S
                for ra, dec, mag in zip(ra_source, dec_source, mag_source):
                    dx_s.append((ra - ra_cat) * 3600)
                    dy_s.append((dec - dec_cat) * 3600)
                    m0_s.append(mag_cat)
                    dm_s.append(mag - mag_cat)

        dx_s = numpy.array(dx_s)
        dy_s = numpy.array(dy_s)
        m0_s = numpy.array(m0_s)
        dm_s = numpy.array(dm_s)

        return dx_m, dy_m, dx_s, dy_s, m0_m, dm_m, m0_s, dm_s

    def makeFluxStat(self, allMat, allSource, calibDic, ffpDic, wcsDic):

        x = list()
        y = list()
        ra = list()
        dec = list()
        id = list()
        for ss in allMat:
            Sxx = 0.0
            Sx  = 0.0
            S   = 0.0
            Sr  = 0.0
            Sd  = 0.0
            if len(ss) > 2:
                for j in range(1,len(ss)):
                    iexp = ss[j].getExp()
                    ichip = ss[j].getChip()
                    if ss[j].getFlux() > 0.0:
                        mag = calibDic[iexp][ichip].getMagnitude(ss[j].getFlux())
                        err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                        xs, ys = ss[j].getX(), ss[j].getY()
                        mcor = ffpDic[iexp][ichip].eval(xs, ys)
                        jcor = -2.5*math.log10(wcsDic[iexp][ichip].pixArea(afwGeom.Point2D(xs, ys)) / wcsDic[iexp][ichip].pixelScale().asDegrees()**2)
                        Sxx += (mag+mcor+jcor)*(mag+mcor+jcor) / (err*err)
                        Sx  += (mag+mcor+jcor) / (err*err)
                        S   += 1. / (err*err)
                        Sr += ss[j].getRa().asDegrees()
                        Sd += ss[j].getDec().asDegrees()
                avg = Sx / S
                sig = math.sqrt(Sxx/S - avg*avg)
                x.append(avg)
                y.append(sig)
                ra.append(Sr / (len(ss)-1))
                dec.append(Sd / (len(ss)-1))

        for ss in allSource:
            Sxx = 0.0
            Sx = 0.0
            S  = 0.0
            Sr  = 0.0
            Sd  = 0.0
            for j in range(1,len(ss)):
                iexp = ss[j].getExp()
                ichip = ss[j].getChip()
                #print iexp, ichip, calibDic[iexp][ichip].getFluxMag0()[0], ss[j].getFlux()
                if calibDic[iexp][ichip].getFluxMag0()[0] > 0 and ss[j].getFlux() > 0.0:
                    mag = calibDic[iexp][ichip].getMagnitude(ss[j].getFlux())
                    err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                    xs, ys = ss[j].getX(), ss[j].getY()
                    mcor = ffpDic[iexp][ichip].eval(xs, ys)
                    jcor = -2.5*math.log10(wcsDic[iexp][ichip].pixArea(afwGeom.Point2D(xs, ys)) / wcsDic[iexp][ichip].pixelScale().asDegrees()**2)
                    Sxx += (mag+mcor+jcor)*(mag+mcor+jcor) / (err*err)
                    Sx  += (mag+mcor+jcor) / (err*err)
                    S   += 1. / (err*err)
                    Sr += ss[j].getRa().asDegrees()
                    Sd += ss[j].getDec().asDegrees()
            if S > 0:
                try:
                    avg = Sx / S
                    sig = math.sqrt(Sxx/S - avg*avg)
                    x.append(avg)
                    y.append(sig)
                    ra.append(Sr / (len(ss)-1))
                    dec.append(Sd / (len(ss)-1))
                except Exception, e:
                    #print Sxx, S, avg, Sxx/S - avg*avg, len(ss)-1
                    pass

        if True:
            plt.clf()
            plt.plot(x, y, ',', markeredgewidth=0)
            plt.xlim(15, 25)
            plt.ylim(0.0, 0.20)
            plt.plot([15, 25], [0.01, 0.01], 'k--')
            plt.xlabel('mag (avg)')
            plt.ylabel('RMS')
            #plt.title('r-band')
            plt.savefig('fluxMean.png')
        else:
            for r, d, m, dm in zip(ra, dec, x, y):
                print '%9.5f %9.5f %7.4f %7.4f' % (r, d, m ,dm)

    def plotPos(self, dx_m, dy_m, dx_s, dy_s):

        x_std_m, x_mean_m, x_n_m = self.clippedStd(dx_m, 3)
        y_std_m, y_mean_m, y_n_m = self.clippedStd(dy_m, 3)
        x_std_s, x_mean_s, x_n_s = self.clippedStd(dx_s, 3)
        y_std_s, y_mean_s, y_n_s = self.clippedStd(dy_s, 3)

        plt.clf()
        plt.rc('text', usetex=True)

        plt.subplot2grid((5,6), (1,0), colspan=4, rowspan=4)
        plt.plot(dx_m, dy_m, 'g,', markeredgecolor='green')
        plt.plot(dx_s, dy_s, 'r,', markeredgecolor='red')
        plt.xlim(-0.5, 0.5)
        plt.ylim(-0.5, 0.5)
        plt.xlabel(r'$\Delta\alpha$ (arcsec)')
        plt.ylabel(r'$\Delta\delta$ (arcsec)')

        bins = numpy.arange(-0.5, 0.5, 0.01) + 0.005

        ax = plt.subplot2grid((5,6),(0,0), colspan=4)
        plt.hist([dx_m, dx_s], bins=bins, normed=False, histtype='step', color=['green', 'red'])
        plt.text(0.75, 0.7, r"$\sigma=$%5.3f" % (x_std_m), transform=ax.transAxes, color='green')
        plt.text(0.75, 0.5, r"$\sigma=$%5.3f" % (x_std_s), transform=ax.transAxes, color='red')
        gauss = mlab.normpdf(bins, x_mean_m, x_std_m)
        plt.plot(bins, gauss*x_n_m*0.01, 'g:')
        gauss = mlab.normpdf(bins, x_mean_s, x_std_s)
        plt.plot(bins, gauss*x_n_s*0.01, 'r:')
        plt.xlim(-0.5, 0.5)

        ax = plt.subplot2grid((5,6),(1,4), rowspan=4)
        plt.hist(dy_m, bins=bins, normed=False, orientation='horizontal', histtype='step', color='green')
        plt.hist(dy_s, bins=bins, normed=False, orientation='horizontal', histtype='step', color='red')
        plt.text(0.7, 0.25, r"$\sigma=$%5.3f" % (y_std_m), rotation=270, transform=ax.transAxes, color='green')
        plt.text(0.5, 0.25, r"$\sigma=$%5.3f" % (y_std_s), rotation=270, transform=ax.transAxes, color='red')
        gauss = mlab.normpdf(bins, y_mean_m, y_std_m)
        plt.plot(gauss*y_n_m*0.01, bins, 'g:')
        gauss = mlab.normpdf(bins, y_mean_s, y_std_s)
        plt.plot(gauss*y_n_s*0.01, bins, 'r:')
        plt.xticks(rotation=270)
        plt.yticks(rotation=270)
        plt.ylim(-0.5, 0.5)

        plt.savefig("posScatter.png")

    def plotFlux(self, m0_m, dm_m, m0_s, dm_s):

        m0Sub = m0_s[m0_s < self.config.maxMag]
        dmSub = dm_s[m0_s < self.config.maxMag]
        
        mag_std_m, mag_mean_m, mag_n_m = self.clippedStd(dm_m, 3)
        mag_std_s, mag_mean_s, mag_n_s = self.clippedStd(dm_s, 3)
        mag_std_sub, mag_mean_sub, mag_n_sub = self.clippedStd(dmSub, 3)

        bins_m = numpy.arange(-0.25, 0.25, 0.025) + 0.0125
        bins_s = numpy.arange(-0.25, 0.25, 0.005) + 0.0025

        plt.clf()
        plt.rc('text', usetex=True)

        plt.subplot2grid((5,6),(1,0), colspan=4, rowspan=4)
        plt.plot(m0_m, dm_m, 'g,', markeredgecolor='green')
        plt.plot(m0_s, dm_s, 'r,', markeredgecolor='red')
        plt.plot(m0Sub, dmSub, 'b,', markeredgecolor='blue')
        plt.plot([15,25], [0,0], 'k--')
        plt.xlim(15, 25)
        plt.ylim(-0.25, 0.25)
        plt.plot([15, 25], [-0.01, -0.01], 'k--')
        plt.plot([15, 25], [+0.01, +0.01], 'k--')
        plt.plot([self.config.maxMag, self.config.maxMag], plt.ylim(), 'b:')
        plt.xlabel(r'$m_{cat}$ (mag)')
        plt.ylabel(r'$\Delta m$ (mag)')

        ax = plt.subplot2grid((5,6),(1,4), rowspan=4)
        plt.hist(dm_s, bins=bins_s, normed=False, orientation='horizontal', histtype='step', color='red')
        plt.hist(dm_m, bins=bins_m, normed=False, orientation='horizontal', histtype='step', color='green')
        plt.hist(dmSub, bins=bins_s, normed=False, orientation='horizontal', histtype='step', color='blue')
        plt.text(0.7, 0.25, r"$\sigma=$%5.3f" % (mag_std_m), rotation=270, transform=ax.transAxes, color='green')
        plt.text(0.5, 0.25, r"$\sigma=$%5.3f" % (mag_std_s), rotation=270, transform=ax.transAxes, color='red')
        plt.text(0.3, 0.25, r"$\sigma=$%5.3f" % (mag_std_sub), rotation=270, transform=ax.transAxes, color='blue')
        gauss = mlab.normpdf(bins_m+0.0125, mag_mean_m, mag_std_m)
        plt.plot(gauss*mag_n_m*0.025, bins_m+0.0125, 'g:')
        gauss = mlab.normpdf(bins_s+0.0025, mag_mean_s, mag_std_s)
        plt.plot(gauss*mag_n_s*0.005, bins_s+0.0025, 'r:')
        gauss = mlab.normpdf(bins_s+0.0025, mag_mean_sub, mag_std_sub)
        plt.plot(gauss*mag_n_sub*0.005, bins_s+0.0025, 'b:')
        plt.xticks(rotation=270)
        plt.yticks(rotation=270)
        plt.ylim(-0.25, 0.25)

        plt.savefig("fluxScatter.png")

    def plotPosAsMag(self, m0_s, dx_s, dy_s):

        plt.subplot2grid((2,5),(0,0), colspan=4)
        plt.plot(m0_s, dx_s, ',', markeredgewidth=0)
        plt.plot([15,25], [0,0], 'k--')
        plt.plot([15,25], [0.01,0.01], 'k--')
        plt.plot([15,25], [-0.01,-0.01], 'k--')
        plt.xlabel('mag')
        plt.ylabel(r'$\Delta\alpha$ (arcsec)')
        plt.xlim(15, 25)
        plt.ylim(-0.15, 0.15)

        mlim = 24
        ax = plt.subplot2grid((2,5),(0,4))
        bins = numpy.arange(-0.15, 0.15, 0.01) + 0.005
        plt.hist(dx_s[m0_s<mlim], bins=bins, normed=False, orientation='horizontal', histtype='step')
        std, mean, n = self.clippedStd(dx_s[m0_s<mlim], 3)
        plt.text(0.7, 0.3, r"$\sigma=$%5.3f" % (std), rotation=270, transform=ax.transAxes, fontsize=10)
        bins = numpy.arange(-0.15, 0.15, 0.001) + 0.0005
        gauss = mlab.normpdf(bins, mean, std)
        plt.plot(gauss*n*0.01, bins)
        plt.ylim(-0.15, 0.15)
        plt.xticks(rotation=270, fontsize=10)
        plt.yticks(rotation=270, fontsize=10)

        plt.subplot2grid((2,5),(1,0), colspan=4)
        plt.plot(m0_s, dy_s, ',', markeredgewidth=0)
        plt.plot([15,25], [0,0], 'k--')
        plt.plot([15,25], [0.01,0.01], 'k--')
        plt.plot([15,25], [-0.01,-0.01], 'k--')
        plt.xlabel('mag')
        plt.ylabel(r'$\Delta\delta$ (arcsec)')
        plt.xlim(15, 25)
        plt.ylim(-0.15, 0.15)

        ax = plt.subplot2grid((2,5),(1,4))
        bins = numpy.arange(-0.15, 0.15, 0.01) + 0.005
        plt.hist(dy_s[m0_s<mlim], bins=bins, normed=False, orientation='horizontal', histtype='step')
        std, mean, n = self.clippedStd(dy_s[m0_s<mlim], 3)
        plt.text(0.7, 0.3, r"$\sigma=$%5.3f" % (std), rotation=270, transform=ax.transAxes, fontsize=10)
        bins = numpy.arange(-0.15, 0.15, 0.001) + 0.0005
        gauss = mlab.normpdf(bins, mean, std)
        plt.plot(gauss*n*0.01, bins)
        plt.ylim(-0.15, 0.15)
        plt.xticks(rotation=270, fontsize=10)
        plt.yticks(rotation=270, fontsize=10)

        plt.savefig('posAsMag.png')

    def writeCatalog(self, allSource, wcsDic, calibDic, ffpDic):
        num = sum(sum(1 for src in ss if src.getExp() >=0 and src.getChip() >= 0) for ss in allSource)

        import pyfits

        schema = pyfits.ColDefs([pyfits.Column(name="id", format="K"),
                                 pyfits.Column(name="ra", format="D"),
                                 pyfits.Column(name="dec", format="D"),
                                 pyfits.Column(name="mag", format="E"),
                                 pyfits.Column(name="err", format="E"),
                                 pyfits.Column(name="corr", format="E"),
                                 ])

        outHdu = pyfits.new_table(schema, nrows=num)
        outData = outHdu.data

        i = 0
        for ss in allSource:
            for src in ss:
                iexp = src.getExp()
                ichip = src.getChip()
                if iexp < 0 or ichip < 0:
                    continue

                outData.id[i] = src.getId()
                x, y = src.getX(), src.getY()
                ra, dec = wcsDic[iexp][ichip].pixelToSky(x, y).getPosition()
                outData.ra[i] = ra
                outData.dec[i] = dec
                fluxMag0 = calibDic[iexp][ichip].getFluxMag0()[0]
                flux = src.getFlux()
                if flux > 0 and fluxMag0 > 0:
                    mcor = ffpDic[iexp][ichip].eval(x, y)
                    jcor = -2.5*math.log10(wcsDic[iexp][ichip].pixArea(afwGeom.Point2D(x, y)) / wcsDic[iexp][ichip].pixelScale().asDegrees()**2)
                    outData.mag[i] = -2.5*math.log10(flux/fluxMag0) + mcor + jcor
                    outData.err[i] = 2.5/math.log(10) * src.getFluxErr() / flux
                    outData.corr[i] = mcor + jcor
                    i += 1

        outHdu.writeto("catalog_check.fits", clobber=True)


    def check(self, dataRefList, ct=None, debug=False, verbose=False):
        ccdSet = self.readCcd(dataRefList)
        self.removeNonExistCcd(dataRefList, ccdSet)

        sourceSet = measMosaic.SourceGroup()
        matchList = measMosaic.SourceMatchGroup()
        astrom = measAstrom.ANetBasicAstrometryTask(self.config.astrom)
        ssVisit = dict()
        mlVisit = dict()
        dataRefListUsed = list()
        wcsDic = dict()
        calibDic = dict()
        ffpDic = dict()
        for dataRef in dataRefList:
            if not ssVisit.has_key(dataRef.dataId['visit']):
                ssVisit[dataRef.dataId['visit']] = list()
                mlVisit[dataRef.dataId['visit']] = list()
                wcsDic[dataRef.dataId['visit']] = dict()
                calibDic[dataRef.dataId['visit']] = dict()
                ffpDic[dataRef.dataId['visit']] = dict()

            try:
                if not dataRef.datasetExists('src'):
                    raise RuntimeError("no data for src %s" % (dataRef.dataId))

                if not dataRef.datasetExists('wcs'):
                    raise RuntimeError("no data for wcs %s" % (dataRef.dataId))

                if not dataRef.datasetExists('fcr'):
                    raise RuntimeError("no data for fcr %s" % (dataRef.dataId))

                md = dataRef.get('wcs_md')
                wcs = afwImage.makeWcs(md)

                md = dataRef.get('calexp_md')
                filterName = afwImage.Filter(md).getName()

                md = dataRef.get('fcr_md')
                ffp = measMosaic.FluxFitParams(md)
                calib = afwImage.Calib(md)

                sources = dataRef.get('src',
                              flags=afwTable.SOURCE_IO_NO_FOOTPRINTS,
                              immediate=True)

                icSrces = dataRef.get('icSrc',
                              flags=afwTable.SOURCE_IO_NO_FOOTPRINTS,
                              immediate=True)
                packedMatches = dataRef.get('icMatch')
                matches = astrom.joinMatchListWithCatalog(packedMatches, icSrces)
                matches[0].first.schema.getAliasMap().set("flux", matches[0].first.schema.join(
                                                          filterName, "flux"))
                matches[0].first.schema.getAliasMap().set("fluxSigma", matches[0].first.schema.join(
                                                          filterName, "fluxSigma"))

                matches = [m for m in matches if m.first != None]
                if ct != None and len(matches) != 0:
                    refSchema = matches[0].first.schema
                    key_p = refSchema.find(refSchema.join(ct.primary, "flux")).key
                    key_s = refSchema.find(refSchema.join(ct.secondary, "flux")).key
                    key_f = refSchema.find("flux").key
                    for m in matches:
                        refFlux1 = m.first.get(key_p)
                        refFlux2 = m.first.get(key_s)
                        refMag1 = -2.5*math.log10(refFlux1)
                        refMag2 = -2.5*math.log10(refFlux2)
                        refMag = ct.transformMags(refMag1, refMag2)
                        refFlux = math.pow(10.0, -0.4*refMag)
                        if refFlux == refFlux:
                            m.first.set(key_f, refFlux)
                        else:
                            m.first = None
                sources = self.selectStars(sources)
                matches = self.selectStars(matches, True)
            except Exception, e:
                print "Failed to read: %s" % (e)
                sources = None
                continue

            if sources != None:
                for s in sources:
                    if numpy.isfinite(s.getRa().asDegrees()): # get rid of NaN
                        src = measMosaic.Source(s)
                        src.setExp(dataRef.dataId['visit'])
                        src.setChip(dataRef.dataId['ccd'])
                        ssVisit[dataRef.dataId['visit']].append(src)
                for m in matches:
                    if m.first != None and m.second != None:
                        match = measMosaic.SourceMatch(measMosaic.Source(m.first, wcs),
                                                       measMosaic.Source(m.second))
                        match.second.setExp(dataRef.dataId['visit'])
                        match.second.setChip(dataRef.dataId['ccd'])
                        mlVisit[dataRef.dataId['visit']].append(match)
                wcsDic[dataRef.dataId['visit']][dataRef.dataId['ccd']] = wcs
                calibDic[dataRef.dataId['visit']][dataRef.dataId['ccd']] = calib
                ffpDic[dataRef.dataId['visit']][dataRef.dataId['ccd']] = ffp
                dataRefListUsed.append(dataRef)

        for visit in ssVisit.keys():
            sourceSet.push_back(ssVisit[visit])
            matchList.push_back(mlVisit[visit])

        d_lim = afwGeom.Angle(self.config.radXMatch, afwGeom.arcseconds)
        nbrightest = self.config.nBrightest

        allMat, allSource = self.mergeCatalog(sourceSet, matchList, ccdSet, d_lim)
 
        dx_m, dy_m, dx_s, dy_s, m0_m, dm_m, m0_s, dm_s  = self.makeDiffPosFlux(allMat, allSource, wcsDic, calibDic, ffpDic)
        self.plotFlux(m0_m, dm_m, m0_s, dm_s)
        self.makeFluxStat(allMat, allSource, calibDic, ffpDic, wcsDic)
        self.plotPos(dx_m, dy_m, dx_s, dy_s)
        self.plotPosAsMag(m0_s, dx_s, dy_s)
        self.writeCatalog(allSource, wcsDic, calibDic, ffpDic)

    def run(self, camera, butler, dataRefList, debug):

        colorterms = ColortermLibraryConfig()
        name = os.path.join(os.environ["OBS_SUBARU_DIR"], "config", camera, "colorterms.py")
        if os.path.exists(name):
            colorterms.load(name)

        filters = list()
        for dataRef in dataRefList:
            if not dataRef.dataId['filter'] in filters:
                filters.append(dataRef.dataId['filter'])

        if len(filters) != 1:
            self.log.fatal("There are %d filters in input frames" % len(filters))
            return None

        ct = colorterms.selectColorTerm(butler.mapper.filters[filters[0]])

        return self.check(dataRefList, ct, debug)

if __name__ == '__main__':

    CheckMosaicTask.parseAndRun()
