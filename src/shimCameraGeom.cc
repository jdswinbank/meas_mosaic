#include "lsst/meas/mosaic/shimCameraGeom.h"

namespace lsst {
namespace meas {
namespace mosaic {

int getNQuarter(CONST_PTR(afw::cameraGeom::Detector) det) {
    return det->getOrientation().getNQuarter() % 4;
}

afw::geom::Angle getYaw(CONST_PTR(afw::cameraGeom::Detector) det) {
    afw::geom::Angle deg = det->getOrientation().getYaw() - getNQuarter(det)*90*afw::geom::degrees;
    deg.wrapCtr();
    return deg;
}

afw::geom::LinearTransform makeScalingMmToPx(afw::geom::Extent2D const pSize) {
    return afw::geom::LinearTransform::makeScaling(1.0/pSize.getX(), 1.0/pSize.getY());
}

afw::geom::Point2D getCenterInFpPixels(CONST_PTR(afw::cameraGeom::Detector) det) {
    auto scaling = makeScalingMmToPx(det->getPixelSize());
    return scaling(det->getCenter(afw::cameraGeom::FOCAL_PLANE).getPoint());
}

afw::geom::Point2D getCenterInDetectorPixels(CONST_PTR(afw::cameraGeom::Detector) det) {
    auto center = det->getCenter(afw::cameraGeom::PIXELS).getPoint();
    if ((getNQuarter(det) % 2) != 0) {
        return afw::geom::Point2D(center.getY(), center.getX());
    } else {
        return center;
    }
}

int getWidth(CONST_PTR(afw::cameraGeom::Detector) det) {
    if ((getNQuarter(det) % 2) != 0) {
        return det->getBBox().getHeight();
    } else {
        return det->getBBox().getWidth();
    }
}

int getHeight(CONST_PTR(afw::cameraGeom::Detector) det) {
    if ((getNQuarter(det) % 2) != 0) {
        return det->getBBox().getWidth();
    } else {
        return det->getBBox().getHeight();
    }
}

afw::geom::Point2D detPxToFpPx(CONST_PTR(afw::cameraGeom::Detector) det, afw::geom::Point2D const detPt) {
    afw::geom::Point2D point;
    if (getNQuarter(det) == 3) {
        point = afw::geom::Point2D(det->getBBox().getWidth() - detPt.getY() - 1, detPt.getX());
    } else if (getNQuarter(det) == 2) {
        point = afw::geom::Point2D(det->getBBox().getWidth() - detPt.getX() - 1,
                                   det->getBBox().getHeight() - detPt.getY() - 1);
    } else if (getNQuarter(det) == 1) {
        point = afw::geom::Point2D(detPt.getY(), det->getBBox().getHeight() - detPt.getX() - 1);
    } else {
        point = detPt;
    }
    auto scaling = makeScalingMmToPx(det->getPixelSize());
    return scaling(det->transform(det->makeCameraPoint(point, afw::cameraGeom::PIXELS),
                                  afw::cameraGeom::FOCAL_PLANE).getPoint());
}

}}} // lsst::meas::mosaic
