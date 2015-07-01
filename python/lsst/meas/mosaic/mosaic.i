// Enable ndarray's NumPy typemaps; types are declared in %included files.
%{
#define PY_ARRAY_UNIQUE_SYMBOL LSST_AFW_MATH_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}
%init %{
    import_array();
%}
%include "ndarray.i"

%{
#include "lsst/meas/mosaic/mosaicfit.h"
#include "lsst/meas/mosaic/fluxfit.h"
#include "lsst/meas/mosaic/spatialCellSource.h"
#include "lsst/meas/mosaic/matches.h"
%}

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"
%include "std_pair.i"

%template(map_int_float) std::map<int, float>;

%import "lsst/afw/cameraGeom/cameraGeomLib.i"
%shared_ptr(lsst::meas::mosaic::Source);
%shared_ptr(lsst::meas::mosaic::Coeff);
%shared_ptr(lsst::meas::mosaic::KDTree);
%shared_ptr(lsst::meas::mosaic::Obs);
%shared_ptr(lsst::meas::mosaic::FluxFitParams);
%shared_ptr(lsst::meas::mosaic::SpatialCellSource);

%inline %{
  PTR(lsst::meas::mosaic::SpatialCellSource) cast_SpatialCellSource(PTR(lsst::afw::math::SpatialCellCandidate) candidate) {
    return boost::dynamic_pointer_cast<lsst::meas::mosaic::SpatialCellSource>(candidate);
  }
%}

%include "lsst/meas/mosaic/mosaicfit.h"
%include "lsst/meas/mosaic/fluxfit.h"
%include "lsst/meas/mosaic/spatialCellSource.h"
%include "lsst/meas/mosaic/matches.h"

%template(SourceSet) std::vector<PTR(lsst::meas::mosaic::Source)>;
%template(SourceGroup) std::vector<std::vector<PTR(lsst::meas::mosaic::Source)> >;
%template(SourceMatch) std::pair<PTR(lsst::meas::mosaic::Source), PTR(lsst::meas::mosaic::Source)>;
%template(SourceMatchSet) std::vector<lsst::meas::mosaic::SourceMatch>;
%template(SourceMatchGroup) std::vector<std::vector<lsst::meas::mosaic::SourceMatch> >;

%template(WcsDic) std::map<int, lsst::afw::image::Wcs::Ptr>;
%template(CcdSet) std::map<int, PTR(lsst::afw::cameraGeom::Detector)>;
%template(CoeffSet) std::map<int, lsst::meas::mosaic::Coeff::Ptr>;
%template(ObsVec) std::vector<lsst::meas::mosaic::Obs::Ptr>;
%template(FfpSet) std::map<int, lsst::meas::mosaic::FluxFitParams::Ptr>;
