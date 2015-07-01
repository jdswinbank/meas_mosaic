#include <stdio.h>
#include <math.h>
#include <vector>
#include <lsst/afw/image/Wcs.h>
#include <lsst/afw/image/Exposure.h>
#include "lsst/meas/mosaic/mosaicfit.h"

using namespace lsst::meas::mosaic;
using namespace lsst::afw::detection;

int main(int argc, char **argv)
{
    char buf[BUFSIZ];

    int order = 9;
    Poly::Ptr p = Poly::Ptr(new Poly(order));
    int nexp = 5;
    int ncoeff = p->ncoeff;
    int nchip = 104;

    CoeffSet coeffs;
    FILE *fp = fopen("coeffs.dat", "rt");
    for (int i = 0; i < nexp; i++) {
	Coeff::Ptr c = Coeff::Ptr(new Coeff(p));
	int k;
	fgets(buf, BUFSIZ, fp);
	sscanf(buf, "%d %lf %lf", &k, &c->A, &c->D);
	for (int j = 0; j < ncoeff; j++) {
	    fgets(buf, BUFSIZ, fp);
	    sscanf(buf, "%d %lf %lf %lf %lf", &k, &c->a[j], &c->b[j], &c->ap[j], &c->bp[j]);
	}
	coeffs.push_back(c);
    }
    fclose(fp);

    double u, v, xi, eta;
    u = 18000.;
    v = 18000.;
    coeffs[0]->uvToXiEta(u, v, &xi, &eta);
    std::cout << u << " " << v << " " << xi << " " << eta << std::endl;
    coeffs[0]->xietaToUV(xi, eta, &u, &v);
    std::cout << u << " " << v << " " << xi << " " << eta << std::endl;

    CcdSet ccdSet;
    fp = fopen("ccd.dat", "rt");
    for (int i = 0; i < nchip; i++) {
	int k;
	double x, y, t;
	fgets(buf, BUFSIZ, fp);
	sscanf(buf, "%d %lf %lf %lf", &k, &x, &y, &t);
	lsst::afw::geom::PointD center = 
	    lsst::afw::geom::Point2D(x, y);
	lsst::afw::cameraGeom::Orientation orientation(0, 0.0, 0.0, t);
	PTR(lsst::afw::cameraGeom::Detector) ccd = PTR(lsst::afw::cameraGeom::Detector)(new lsst::afw::cameraGeom::Ccd(k, 0.168));
	ccd->setCenter(lsst::afw::cameraGeom::FpPoint(center));
	ccd->setOrientation(orientation);
	ccdSet.push_back(ccd);
    }
    fclose(fp);

    for (size_t j = 0; j < 1; j++) {
	for (size_t i = 99; i < 100; i++) {
	    //    for (size_t j = 0; j < coeffs.size(); j++) {
	    //	for (size_t i = 0; i < ccdSet.size(); i++) {
	    Coeff::Ptr c = convertCoeff(coeffs[j], ccdSet[i]);
	    lsst::afw::image::TanWcs::ConstPtr wcs = wcsFromCoeff(c);
	    lsst::afw::image::Exposure<int> exposure(lsst::afw::geom::Extent2I(0,0), wcs);

	    lsst::afw::coord::Coord::Ptr cp = wcs->pixelToSky(2000, 4000);
	    std::cout << cp->getLongitude().asDegrees() << " " << cp->getLatitude().asDegrees() << std::endl;
	    lsst::afw::geom::PointD pt = wcs->skyToPixel(cp);
	    std::cout << pt[0] << " " << pt[1] << std::endl;

	    char fname[BUFSIZ];
	    sprintf(fname, "wcs%ld%03ld.fits", j, i);
	    exposure.writeFits(fname);
	}
    }

    return 0;
}
