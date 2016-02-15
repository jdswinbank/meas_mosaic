// -*- lsst-c++ -*-
%define fitLib_DOCSTRING
"
Python interface to lsst::meas::mosaic
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.mosaic", docstring=fitLib_DOCSTRING) mosaicLib

%{
#include "lsst/afw/image.h"
#include "lsst/pex/logging.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/math.h"
#include "lsst/afw/table.h"
#include "lsst/afw/geom/polygon/Polygon.h"
%}

%include "lsst/p_lsstSwig.i"

%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/math/mathLib.i"
%import "lsst/afw/table/tableLib.i"

%include "mosaic.i"
