// -*- lsst-c++ -*-
#if !defined(HSC_MEAS_MOSAIC_H)
#define HSC_MEAS_MOSAIC_H

#include <vector>
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/table.h"
#include "lsst/utils/ieee.h"
#include "boost/enable_shared_from_this.hpp"

namespace lsst {
    namespace meas {
	namespace mosaic {

            class Source {
            public:
                enum { UNSET = -1 };
                typedef long IdType;
                typedef int ChipType;
                typedef int ExpType;
                explicit Source(lsst::afw::table::SourceRecord const& record) :
                    _id(record.getId()), _chip(UNSET), _exp(UNSET), _sky(record.getRa(), record.getDec()),
                    _pixels(record.getX(), record.getY()), _flux(record.getCalibFlux()), _err(record.getCalibFluxErr()),
		    _xerr(sqrt(record.getCentroidErr()(0,0))), _yerr(sqrt(record.getCentroidErr()(1,1))),
                    _astromBad(record.getCentroidFlag() | record.getCalibFluxFlag()) {}
                Source(lsst::afw::table::SimpleRecord const& record, lsst::afw::image::Wcs const& wcs) :
                    _id(record.getId()), _chip(UNSET), _exp(UNSET), _sky(record.getRa(), record.getDec()),
                    _pixels(wcs.skyToPixel(_sky)),
                    _flux(record.get(record.getSchema().find<double>("flux").key)),
                    _err(std::numeric_limits<double>::quiet_NaN()),
                    _xerr(std::numeric_limits<double>::quiet_NaN()),
                    _yerr(std::numeric_limits<double>::quiet_NaN()),
                    _astromBad(!lsst::utils::isfinite(_flux))
                    {
                        try {
                            _err = record.get(record.getSchema().find<double>("fluxSigma").key);
                        } catch (pex::exceptions::NotFoundError const&) {
			    // flux.err is not availabe. sqrt of flux will be used
			    // 1.E-08 is a scaling factor to make mag_err ~ 0.1mag for 15 mag
			    // This is a good approximation for 2MASS J-band
			    _err = sqrt(record.get(record.getSchema().find<double>("flux").key)*1.E-08);
			}
                    }
                Source(lsst::afw::coord::Coord coord, double flux=std::numeric_limits<double>::quiet_NaN()) :
                    _id(-1), _chip(UNSET), _exp(UNSET), _sky(coord),
                    _pixels(lsst::afw::geom::Point2D(std::numeric_limits<double>::quiet_NaN(),
                                                     std::numeric_limits<double>::quiet_NaN())),
                    _flux(flux), _err(std::numeric_limits<double>::quiet_NaN()),
                    _xerr(std::numeric_limits<double>::quiet_NaN()),
                    _yerr(std::numeric_limits<double>::quiet_NaN()),
                    _astromBad(false) {}

		// This ctor is for python pickling
		Source(IdType id, ChipType chip, ExpType exp,
		    double ra, double dec,
		    double x, double xerr, double y, double yerr,
		    double flux, double fluxerr,
		    bool astromBad
		):
		    _id(id), _chip(chip), _exp(exp),
		    _sky(lsst::afw::geom::Angle(ra, lsst::afw::geom::degrees), lsst::afw::geom::Angle(dec, lsst::afw::geom::degrees)),
		    _pixels(x, y),
		    _flux(flux), _err(fluxerr),
		    _xerr(xerr), _yerr(yerr),
		    _astromBad(astromBad)
		{}

                IdType getId() const { return _id; }
                ChipType getChip() const { return _chip; }
                ExpType getExp() const { return _exp; }
                lsst::afw::coord::Coord getSky() const { return _sky; }
                lsst::afw::geom::Angle getRa() const { return getSky().getLongitude(); }
                lsst::afw::geom::Angle getDec() const { return getSky().getLatitude(); }
                lsst::afw::geom::Point2D getPixels() const { return _pixels; }
                double getX() const { return getPixels().getX(); }
                double getY() const { return getPixels().getY(); }
		double getXErr() const { return _xerr; }
		double getYErr() const { return _yerr; }
                double getFlux() const { return _flux; }
                double getFluxErr() const { return _err; }
                bool getAstromBad() const { return _astromBad; }

                void setChip(ChipType chip) { _chip = chip; }
                void setExp(ExpType exp) { _exp = exp; }
                void setFlux(double value) { _flux = value; }

            private:
                IdType _id;                       // Identifier
                ChipType _chip;                   // Chip identifier
                ExpType _exp;                     // Exposure identifier
                lsst::afw::coord::Coord _sky;     // Sky coordinates
                lsst::afw::geom::Point2D _pixels; // Pixel coordinates
                double _flux;                     // Flux
                double _err;			  // Flux Err
		double _xerr;			  // x coordinate error
		double _yerr;			  // y coordinate error
                bool _astromBad;                  // Astrometry bad?
            };

            typedef std::pair<PTR(Source), PTR(Source)> SourceMatch;

            typedef std::vector<PTR(Source)> SourceSet;
	    typedef std::vector<std::vector<PTR(Source)> > SourceGroup;
            typedef std::vector<SourceMatch> SourceMatchSet;
	    typedef std::vector<std::vector<SourceMatch> > SourceMatchGroup;

	    typedef std::map<int, lsst::afw::image::Wcs::Ptr> WcsDic;

	    class Poly {
	    public:
	        typedef boost::shared_ptr<Poly> Ptr;

		int order;
		int ncoeff;
		int *xorder;
		int *yorder;

		Poly(int order);
		~Poly(void);
		Poly(const Poly &p);
		int getIndex(int i, int j);
		int getXorder(int i) { return xorder[i]; }
		int getYorder(int i) { return yorder[i]; }
	    };

	    class Coeff {
	    public:
		typedef boost::shared_ptr<Coeff> Ptr;

		Poly::Ptr p;
		int iexp;
		double *a;
		double *b;
		double *ap;
		double *bp;
		double A;
		double D;
		double x0;
		double y0;

		Coeff(int order);
		Coeff(Poly::Ptr const & p);
		~Coeff(void);
		Coeff(const Coeff &c);
		void show(void);
		void uvToXiEta(double u, double v, double *xi, double *eta);
		void xietaToUV(double xi, double eta, double *u, double *v);
		double get_a(int i) { return a[i]; }
		double get_b(int i) { return b[i]; }
		double get_ap(int i) { return ap[i]; }
		double get_bp(int i) { return bp[i]; }
		void set_a(int i, double v) { a[i] = v; }
		void set_b(int i, double v) { b[i] = v; }
		void set_ap(int i, double v) { ap[i] = v; }
		void set_bp(int i, double v) { bp[i] = v; }
		double xi(double u, double v);
		double eta(double u, double v);
		double dxidu(double u, double v);
		double dxidv(double u, double v);
		double detadu(double u, double v);
		double detadv(double u, double v);
		double detJ(double u, double v);
		int getNcoeff() { return p->ncoeff; }
		double pixelScale(void);

                // below lines for coeffFromTanWcs()
		void set_D(double v) { D = v; }
		void set_A(double v) { A = v; }
		void set_x0(double v) { x0 = v; }
		void set_y0(double v) { y0 = v; }
		void set_iexp(int v) { iexp = v; }
		double get_D() { return D; }
		double get_A() { return A; }
		double get_x0() { return x0; }
		double get_y0() { return y0; }
		int get_iexp() { return iexp; }
	    };

	    class Obs {
	    public:
		typedef boost::shared_ptr<Obs> Ptr;

		double ra, dec;
		double xi, eta;
		double xi_a, xi_d, eta_a, eta_d;
		double xi_A, xi_D, eta_A, eta_D;

		double x, y;
		double u, v;
		double u0, v0;

		double U, V;

		double xi_fit, eta_fit;
		double u_fit, v_fit;

		long id;
		int istar;
		int jstar;	/* index for fit */
		int iexp;	/* exposure id or visit number */
		int ichip;	/* ccdId */
		int jexp;	/* sequential index for exposure used for matrix*/
		int jchip;	/* sequential index for ccd used for matrix*/
		bool good;

		double xerr;
		double yerr;

		double mag;
		double mag0;
		double err;
		double mag_cat;
		double err_cat;

		Obs(int id, double ra, double dec, double x, double y, int ichip, int iexp);
		Obs(int id, double ra, double dec, int ichip, int iexp);
		void setUV(PTR(lsst::afw::cameraGeom::Detector) &ccd, double x0=0, double y0=0);
		void setXiEta(double ra_c, double dec_c);
		void setFitVal(Coeff::Ptr& c, Poly::Ptr p);
		void setFitVal2(Coeff::Ptr& c, Poly::Ptr p);
	    };

	    class KDTree
#if !defined(SWIG)
	      : public boost::enable_shared_from_this<KDTree>
#endif
	    {
	    public:
		typedef boost::shared_ptr<KDTree> Ptr;
		typedef boost::shared_ptr<const KDTree> ConstPtr;

		int depth;
		int axis;
                lsst::afw::geom::Angle location[2];
		lsst::afw::coord::Coord c;
		KDTree::Ptr left;
		KDTree::Ptr right;
                SourceSet set;

                KDTree(SourceSet& s, int depth) { _initializeSources(s, depth); }
		KDTree(PTR(Source) s, int depth);

                KDTree(SourceMatchSet m, int depth) { _initializeMatches(m, depth); }
                KDTree(SourceMatch const& m, int depth);

		~KDTree();
                ConstPtr search(lsst::afw::coord::Coord const& sky) const;
		ConstPtr findSource(Source const& s) const;
                void add(SourceMatch const& m);
		void add(PTR(Source) s,
                         lsst::afw::geom::Angle d_lim=lsst::afw::geom::Angle(0, lsst::afw::geom::degrees));
		int count(void);
		SourceGroup mergeMat() const;
		SourceGroup mergeSource(unsigned int minNumMatch = 2);
		void printMat() const;
		void printSource() const;
		bool isLeaf(void) const { return left == NULL && right == NULL; }
		KDTree::Ptr findNearest(Source const& s);

		double distance(Source const& s) const {
                    return c.angularSeparation(lsst::afw::coord::Coord(s.getRa(), s.getDec())).asDegrees();
                }

            private:
                void _initializeSources(SourceSet& s, int depth);
                void _initializeMatches(SourceMatchSet& m, int depth);
            };

	    typedef std::map<int, PTR(lsst::afw::cameraGeom::Detector)> CcdSet;
	    typedef std::map<int, Coeff::Ptr> CoeffSet;
	    typedef std::vector<Obs::Ptr> ObsVec;

	    int flagSuspect(SourceGroup &allMat,
			    SourceGroup &allSource,
			    WcsDic &wcsDic);

	    KDTree::Ptr kdtreeMat(SourceMatchGroup &matchList);
	    KDTree::Ptr kdtreeSource(SourceGroup const &sourceSet,
				     KDTree::Ptr rootMat,
				     CcdSet &ccdSet,
				     lsst::afw::geom::Angle d_lim);

	    ObsVec obsVecFromSourceGroup(SourceGroup const &all,
					 WcsDic &wcsDic,
					 CcdSet &ccdSet);

	    CoeffSet solveMosaic_CCD_shot(int order,
					  int nmatch,
					  ObsVec &matchVec,
					  WcsDic &wcsDic,
					  CcdSet &ccdSet,
					  bool solveCcd = true,
					  bool allowRotation = true,
					  bool verbose = false,
					  double catRMS = 0.0,
                                          bool writeSnapshots = false,
                                          std::string const & snapshotDir = ".");

	    CoeffSet solveMosaic_CCD(int order,
				     int nmatch,
				     int nsource,
				     ObsVec &matchVec,
				     ObsVec &sourceVec,
				     WcsDic &wcsDic,
				     CcdSet &ccdSet,
				     bool solveCcd = true,
				     bool allowRotation = true,
				     bool verbose = false,
				     double catRMS = 0.0,
                                     bool writeSnapshots = false,
                                     std::string const & snapshotDir = ".");

	    Coeff::Ptr convertCoeff(Coeff::Ptr& coeff,
				    PTR(lsst::afw::cameraGeom::Detector)& ccd);

	    lsst::afw::image::TanWcs::Ptr wcsFromCoeff(Coeff::Ptr& coeff);

            Coeff::Ptr coeffFromTanWcs(lsst::afw::image::Wcs::Ptr& wcs);

	    lsst::afw::image::Image<float>::Ptr
	      getJImg(Coeff::Ptr& coeff,
		      PTR(lsst::afw::cameraGeom::Detector)& ccd);

	    lsst::afw::image::Image<float>::Ptr
	      getJImg(lsst::afw::image::Wcs::Ptr& wcs,
		      int width, int height);

	    lsst::afw::image::Image<float>::Ptr
	      getJImg(lsst::afw::image::Wcs::Ptr& wcs,
		      PTR(lsst::afw::cameraGeom::Detector)& ccd);

            //{
            // Calculate Jacobian correction
            //
            // This needs to be included along with the flux correction
            // when correcting the flux of objects. It corrects for the
            // dimming of the flat-field due to optical scale variations
            // over the field of view.
            double calculateJacobian(
                afw::image::Wcs const& wcs, ///< Astrometric solution
                afw::geom::Point2D const& point ///< Position for correction
                );
            ndarray::Array<double, 1> calculateJacobian(
                afw::image::Wcs const& wcs, ///< Astrometric solution
                std::vector<afw::geom::Point2D> const& points ///< Positions for correction
                );
            ndarray::Array<double, 1> calculateJacobian(
                afw::image::Wcs const& wcs, ///< Astrometric solution
                ndarray::Array<double const, 1> const& x, ///< x positions for correction
                ndarray::Array<double const, 1> const& y  ///< y positions for correction
                );
            //}
    }
  }
}

#endif
