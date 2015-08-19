#include <ctime>
#include <strings.h>
#include "fitsio.h"

#include "lsst/utils/ieee.h"
#include "lsst/meas/mosaic/mosaicfit.h"
#include "lsst/meas/mosaic/snapshot.h"
#include "lsst/meas/mosaic/shimCameraGeom.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/table/Match.h"
#include "boost/make_shared.hpp"
#include "boost/format.hpp"
#include "boost/filesystem/path.hpp"

#include "dynamic_lapack.h"

#define D2R (M_PI/180.)
#define R2D (180./M_PI)

using namespace lsst::meas::mosaic;

#include "Eigen/Core"
#include "Eigen/LU"

Eigen::VectorXd solveMatrix(long size, Eigen::MatrixXd &a_data, Eigen::VectorXd &b_data);

static void decodeSipHeader(CONST_PTR(lsst::daf::base::PropertySet) const& fitsMetadata,
                            std::string const& which,
                            Eigen::MatrixXd *m);


double calXi(double a, double d, double A, double D);
double calXi_a(double a, double d, double A, double D);
double calXi_d(double a, double d, double A, double D);
double calXi_A(double a, double d, double A, double D);
double calXi_D(double a, double d, double A, double D);
double calEta(double a, double d, double A, double D);
double calEta_a(double a, double d, double A, double D);
double calEta_d(double a, double d, double A, double D);
double calEta_A(double a, double d, double A, double D);
double calEta_D(double a, double d, double A, double D);

Poly::Poly(int order) {
    this->order = order;
    this->ncoeff = (order+1) * (order+2) / 2 - 1;
    xorder = new int[ncoeff];
    yorder = new int[ncoeff];

    int k = 0;
    for (int j = 1; j <= order; j++) {
	for (int i = 0; i <= j; i++) {
	    xorder[k] = j - i;
	    yorder[k] = i;
	    k++;
	}
    }
}

Poly::~Poly(void) {
    delete [] xorder;
    delete [] yorder;
}

Poly::Poly(const Poly &p) {
    this->order = p.order;
    this->ncoeff = p.ncoeff;
    xorder = new int[this->ncoeff];
    yorder = new int[this->ncoeff];
    for (int i = 0; i < this->ncoeff; i++) {
	this->xorder[i] = p.xorder[i];
	this->yorder[i] = p.yorder[i];
    }
}

int Poly::getIndex(int i, int j) {
    for (int k = 0; k < this->ncoeff; k++) {
	if (xorder[k] == i &&
	    yorder[k] == j) {
	    return k;
	}
    }

    return -1;
}

Coeff::Coeff(int order) {
    this->p = Poly::Ptr(new Poly(order));
    this->a  = new double[this->p->ncoeff];
    this->b  = new double[this->p->ncoeff];
    this->ap = new double[this->p->ncoeff];
    this->bp = new double[this->p->ncoeff];
    memset(this->a,  0x0, sizeof(double)*this->p->ncoeff);
    memset(this->b,  0x0, sizeof(double)*this->p->ncoeff);
    memset(this->ap, 0x0, sizeof(double)*this->p->ncoeff);
    memset(this->bp, 0x0, sizeof(double)*this->p->ncoeff);
}

Coeff::Coeff(Poly::Ptr const &p) {
    this->p = p;
    this->a  = new double[this->p->ncoeff];
    this->b  = new double[this->p->ncoeff];
    this->ap = new double[this->p->ncoeff];
    this->bp = new double[this->p->ncoeff];
    memset(this->a,  0x0, sizeof(double)*this->p->ncoeff);
    memset(this->b,  0x0, sizeof(double)*this->p->ncoeff);
    memset(this->ap, 0x0, sizeof(double)*this->p->ncoeff);
    memset(this->bp, 0x0, sizeof(double)*this->p->ncoeff);
}

Coeff::~Coeff(void) {
    delete [] this->a;
    delete [] this->b;
    delete [] this->ap;
    delete [] this->bp;
}

Coeff::Coeff(const Coeff &c) {
    this->p = c.p;
    this->a  = new double[this->p->ncoeff];
    this->b  = new double[this->p->ncoeff];
    this->ap = new double[this->p->ncoeff];
    this->bp = new double[this->p->ncoeff];
    for (int i = 0; i < this->p->ncoeff; i++) {
	this->a[i]  = c.a[i];
	this->b[i]  = c.b[i];
	this->ap[i] = c.ap[i];
	this->bp[i] = c.bp[i];
    }
    this->A = c.A;
    this->D = c.D;
    this->x0 = c.x0;
    this->y0 = c.y0;
}

void Coeff::show(void) {
    printf("%12.5e %12.5e\n", this->A, this->D);
    for (int k = 0; k < this->p->ncoeff; k++) {
	printf("%12.5e %12.5e %12.5e %12.5e\n",
	       this->a[k], this->b[k],
	       this->ap[k], this->bp[k]);
    }
}

void Coeff::uvToXiEta(double u, double v, double *xi, double *eta) {
    *xi  = 0.0;
    *eta = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	*xi  += this->a[i] * pow(u, p->xorder[i]) * pow(v, p->yorder[i]);
	*eta += this->b[i] * pow(u, p->xorder[i]) * pow(v, p->yorder[i]);
    }
}

void Coeff::xietaToUV(double xi, double eta, double *u, double *v) {
    Eigen::Matrix2d cd;
    cd << this->a[0], this->a[1], this->b[0], this->b[1];
    double det = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    double U = ( xi * cd(1,1) - eta * cd(0,1)) / det;
    double V = (-xi * cd(1,0) + eta * cd(1,1)) / det;
    *u = U;
    *v = V;
    for (int i = 0; i < this->p->ncoeff; i++) {
	*u += this->ap[i] * pow(U, p->xorder[i]) * pow(V, p->yorder[i]);
	*v += this->bp[i] * pow(U, p->xorder[i]) * pow(V, p->yorder[i]);
    }
}

double Coeff::xi(double u, double v) {
    double xi = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	xi += this->a[i] * pow(u, this->p->xorder[i]) * pow(v, this->p->yorder[i]);
    }
    return xi;
}

double Coeff::eta(double u, double v) {
    double eta = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	eta += this->b[i] * pow(u, this->p->xorder[i]) * pow(v, this->p->yorder[i]);
    }
    return eta;
}

double Coeff::dxidu(double u, double v) {
    double dxi = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	if (this->p->xorder[i]-1 >= 0) {
	    dxi += this->a[i] * this->p->xorder[i] * pow(u, this->p->xorder[i]-1) * pow(v, this->p->yorder[i]);
	}
    }
    return dxi;
}

double Coeff::dxidv(double u, double v) {
    double dxi = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	if (this->p->yorder[i]-1 >= 0) {
	    dxi += this->a[i] * pow(u, this->p->xorder[i]) * this->p->yorder[i] * pow(v, this->p->yorder[i]-1);
	}
    }
    return dxi;
}

double Coeff::detadu(double u, double v) {
    double deta = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	if (this->p->xorder[i]-1 >= 0) {
	    deta += this->b[i] * this->p->xorder[i] * pow(u, this->p->xorder[i]-1) * pow(v, this->p->yorder[i]);
	}
    }
    return deta;
}

double Coeff::detadv(double u, double v) {
    double deta = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	if (this->p->yorder[i]-1 >= 0) {
	    deta += this->b[i] * pow(u, this->p->xorder[i]) * this->p->yorder[i] * pow(v, this->p->yorder[i]-1);
	}
    }
    return deta;
}

double Coeff::detJ(double u, double v) {
    double a = this->dxidu(u, v);
    double b = this->dxidv(u, v);
    double c = this->detadu(u, v);
    double d = this->detadv(u, v);

    return fabs(a*d-b*c);
}

double Coeff::pixelScale(void) {
    return sqrt(fabs(a[0] * b[1] - a[1] * b[0]));
}

Obs::Obs(int id_, double ra_, double dec_, double x_, double y_, int ichip_, int iexp_) :
    ra(ra_),
    dec(dec_),
    xi(std::numeric_limits<double>::quiet_NaN()),
    eta(std::numeric_limits<double>::quiet_NaN()),
    xi_a(std::numeric_limits<double>::quiet_NaN()),
    xi_d(std::numeric_limits<double>::quiet_NaN()),
    eta_a(std::numeric_limits<double>::quiet_NaN()),
    eta_d(std::numeric_limits<double>::quiet_NaN()),
    xi_A(std::numeric_limits<double>::quiet_NaN()),
    xi_D(std::numeric_limits<double>::quiet_NaN()),
    eta_A(std::numeric_limits<double>::quiet_NaN()),
    eta_D(std::numeric_limits<double>::quiet_NaN()),
    x(x_),
    y(y_),
    u(std::numeric_limits<double>::quiet_NaN()),
    v(std::numeric_limits<double>::quiet_NaN()),
    u0(std::numeric_limits<double>::quiet_NaN()),
    v0(std::numeric_limits<double>::quiet_NaN()),
    U(std::numeric_limits<double>::quiet_NaN()),
    V(std::numeric_limits<double>::quiet_NaN()),
    xi_fit(std::numeric_limits<double>::quiet_NaN()),
    eta_fit(std::numeric_limits<double>::quiet_NaN()),
    u_fit(std::numeric_limits<double>::quiet_NaN()),
    v_fit(std::numeric_limits<double>::quiet_NaN()),
    id(id_),
    istar(-1),
    jstar(-2),
    iexp(iexp_),
    ichip(ichip_),
    jexp(-1),
    good(true),
    mag(std::numeric_limits<double>::quiet_NaN()),
    mag0(std::numeric_limits<double>::quiet_NaN()),
    mag_cat(std::numeric_limits<double>::quiet_NaN())
{}


Obs::Obs(int id_, double ra_, double dec_, int ichip_, int iexp_) :
    ra(ra_),
    dec(dec_),
    xi(std::numeric_limits<double>::quiet_NaN()),
    eta(std::numeric_limits<double>::quiet_NaN()),
    xi_a(std::numeric_limits<double>::quiet_NaN()),
    xi_d(std::numeric_limits<double>::quiet_NaN()),
    eta_a(std::numeric_limits<double>::quiet_NaN()),
    eta_d(std::numeric_limits<double>::quiet_NaN()),
    xi_A(std::numeric_limits<double>::quiet_NaN()),
    xi_D(std::numeric_limits<double>::quiet_NaN()),
    eta_A(std::numeric_limits<double>::quiet_NaN()),
    eta_D(std::numeric_limits<double>::quiet_NaN()),
    x(std::numeric_limits<double>::quiet_NaN()),
    y(std::numeric_limits<double>::quiet_NaN()),
    u(std::numeric_limits<double>::quiet_NaN()),
    v(std::numeric_limits<double>::quiet_NaN()),
    u0(std::numeric_limits<double>::quiet_NaN()),
    v0(std::numeric_limits<double>::quiet_NaN()),
    U(std::numeric_limits<double>::quiet_NaN()),
    V(std::numeric_limits<double>::quiet_NaN()),
    xi_fit(std::numeric_limits<double>::quiet_NaN()),
    eta_fit(std::numeric_limits<double>::quiet_NaN()),
    u_fit(std::numeric_limits<double>::quiet_NaN()),
    v_fit(std::numeric_limits<double>::quiet_NaN()),
    id(id_),
    istar(-1),
    jstar(-2),
    iexp(iexp_),
    ichip(ichip_),
    jexp(-1),
    good(true),
    mag(std::numeric_limits<double>::quiet_NaN()),
    mag0(std::numeric_limits<double>::quiet_NaN()),
    mag_cat(std::numeric_limits<double>::quiet_NaN())
{}

void Obs::setUV(PTR(lsst::afw::cameraGeom::Detector) &ccd, double x0, double y0) {
    double cosYaw = std::cos(getYaw(ccd));
    double sinYaw = std::sin(getYaw(ccd));

    this->u0 = this->x * cosYaw - this->y * sinYaw;
    this->v0 = this->x * sinYaw + this->y * cosYaw;

    afw::geom::Point2D xy(this->x, this->y);
    afw::geom::Point2D uv = detPxToFpPx(ccd, xy);

    this->u = uv.getX() + x0;
    this->v = uv.getY() + y0;
}

void Obs::setXiEta(double ra_c, double dec_c) {
    this->xi    = calXi   (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->eta   = calEta  (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->xi_a  = calXi_a (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->xi_d  = calXi_d (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->eta_a = calEta_a(this->ra, this->dec, ra_c, dec_c) * R2D;
    this->eta_d = calEta_d(this->ra, this->dec, ra_c, dec_c) * R2D;
    this->xi_A  = calXi_A (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->xi_D  = calXi_D (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->eta_A = calEta_A(this->ra, this->dec, ra_c, dec_c) * R2D;
    this->eta_D = calEta_D(this->ra, this->dec, ra_c, dec_c) * R2D;
}

void Obs::setFitVal(Coeff::Ptr& c, Poly::Ptr p) {
    this->xi_fit  = 0.0;
    this->eta_fit = 0.0;
    for (int k = 0; k < c->p->ncoeff; k++) {
	this->xi_fit  += c->a[k] * pow(this->u, p->xorder[k]) * pow(this->v, p->yorder[k]);
	this->eta_fit += c->b[k] * pow(this->u, p->xorder[k]) * pow(this->v, p->yorder[k]);
    }
}

void Obs::setFitVal2(Coeff::Ptr& c, Poly::Ptr p) {
    Eigen::Matrix2d cd;
    cd << c->a[0], c->a[1], c->b[0], c->b[1];
    double det = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    double U = ( this->xi * cd(1,1) - this->eta * cd(0,1)) / det;
    double V = (-this->xi * cd(1,0) + this->eta * cd(0,0)) / det;
    this->u_fit = U;
    this->v_fit = V;
    for (int i = 0; i < c->p->ncoeff; i++) {
	this->u_fit += c->ap[i] * pow(U, p->xorder[i]) * pow(V, p->yorder[i]);
	this->v_fit += c->bp[i] * pow(U, p->xorder[i]) * pow(V, p->yorder[i]);
    }
}

struct SourceMatchCmpRa
{
    template<class MatchT>
    bool operator()(MatchT const& lhs, MatchT const& rhs) const
    {
        return lhs.first->getRa() < rhs.first->getRa();
    }
};

class SourceMatchCmpDec
{
public:
    template<class MatchT>
    bool operator()(MatchT const& lhs, MatchT const& rhs) const
    {
        return lhs.first->getDec() < rhs.first->getDec();
    }
};

class SourceCmpRa
{
public:
    template<class SourceT>
    bool operator()(SourceT const& lhs, SourceT const& rhs) const
    {
        return lhs->getRa() < rhs->getRa();
    }
};

class SourceCmpDec
{
public:
    template<class SourceT>
    bool operator()(SourceT const& lhs, SourceT const& rhs ) const
    {
        return lhs->getDec() < rhs->getDec();
    }
};

KDTree::KDTree(PTR(Source) s, int depth) {
    SourceSet set;
    set.push_back(s);
    _initializeSources(set, depth);
}

KDTree::KDTree(SourceMatch const& m, int depth) {
    std::vector<SourceMatch> matches;
    matches.push_back(m);
    _initializeMatches(matches, depth);
}

void KDTree::_initializeSources(SourceSet& s, int depth)
{
    this->depth = depth;
    this->axis = depth % 2;

    if (s.size() == 1) {

	this->location[0] = s[0]->getRa();
	this->location[1] = s[0]->getDec();
	this->c = lsst::afw::coord::Coord(this->location[0], this->location[1]);
	this->set.push_back(s[0]);

	this->left  = KDTree::Ptr();
	this->right = KDTree::Ptr();

    } else {

	if (this->axis == 0)
	    std::sort(s.begin(), s.end(), SourceCmpRa());
	else
	    std::sort(s.begin(), s.end(), SourceCmpDec());

	this->location[0] = s[s.size()/2]->getRa();
	this->location[1] = s[s.size()/2]->getDec();
	this->c = lsst::afw::coord::Coord(this->location[0], this->location[1]);

	this->set.push_back(s[s.size()/2]);

        size_t middle = s.size() / 2;
        SourceSet s_left(s.begin(), s.begin() + middle);
        SourceSet s_right(s.begin() + middle + 1, s.end());

	if (s_left.size() > 0) {
	    this->left  = KDTree::Ptr(new KDTree(s_left,  depth+1));
	} else {
	    this->left  = KDTree::Ptr();
	}

	if (s_right.size() > 0) {
	    this->right = KDTree::Ptr(new KDTree(s_right, depth+1));
	} else {
	    this->right = KDTree::Ptr();
	}

    }
}


void KDTree::_initializeMatches(SourceMatchSet &m, int depth) {
    this->depth = depth;
    this->axis = depth % 2;

    if (m.size() == 1) {

	this->location[0] = m[0].first->getRa();
	this->location[1] = m[0].first->getDec();
	this->c = lsst::afw::coord::Coord(this->location[0], this->location[1]);

	this->set.push_back(m[0].first);
	this->set.push_back(m[0].second);

	this->left  = KDTree::Ptr();
	this->right = KDTree::Ptr();

    } else {
        if (this->axis == 0) {
            std::sort(m.begin(), m.end(), SourceMatchCmpRa());
        } else {
            std::sort(m.begin(), m.end(), SourceMatchCmpDec());
        }

        size_t middle = m.size() / 2;

	this->location[0] = m[middle].first->getRa();
	this->location[1] = m[middle].first->getDec();
	this->c = lsst::afw::coord::Coord(this->location[0], this->location[1]);

	this->set.push_back(m[middle].first);
	this->set.push_back(m[middle].second);

        std::vector<SourceMatch> m_left(m.begin(), m.begin() + middle);
        std::vector<SourceMatch> m_right(m.begin() + middle + 1, m.end());

	if (m_left.size() > 0) {
	    this->left  = KDTree::Ptr(new KDTree(m_left, depth+1));
	} else {
	    this->left  = KDTree::Ptr();
	}

	if (m_right.size() > 0) {
	    this->right = KDTree::Ptr(new KDTree(m_right, depth+1));
	} else {
	    this->right = KDTree::Ptr();
	}

    }
}

KDTree::~KDTree() {
    if (this->left != NULL) {
	this->left.reset();
    }

    if (this->right != NULL) {
	this->right.reset();
    }
}

KDTree::ConstPtr KDTree::search(lsst::afw::coord::Coord const& sky) const {

    lsst::afw::geom::Angle ra  = sky.getLongitude();
    lsst::afw::geom::Angle dec = sky.getLatitude();

    lsst::afw::geom::Angle val;
    if (this->axis == 0)
	val = ra;
    else
	val = dec;

    if (this->set[0]->getRa()  == ra &&
	this->set[0]->getDec() == dec) {
	return shared_from_this();
    } else {
	if (val < this->location[this->axis]) {
	    if (this->left != NULL) {
		return this->left->search(sky);
	    } else {
		return KDTree::Ptr();
	    }
	} else {
	    if (this->right != NULL) {
		return this->right->search(sky);
	    } else {
		return KDTree::Ptr();
	    }
	}
    }
}

void KDTree::add(SourceMatch const& m) {

    lsst::afw::geom::Angle ra  = m.first->getRa();
    lsst::afw::geom::Angle dec = m.first->getDec();

    lsst::afw::geom::Angle val;
    if (this->axis == 0)
	val = ra;
    else
	val = dec;

    if (this->set[0]->getRa()  == ra &&
	this->set[0]->getDec() == dec) {
	this->set.push_back(m.second);
    } else {
	if (val < this->location[this->axis]) {
	    if (this->left != NULL) {
		this->left->add(m);
	    } else {
		this->left = KDTree::Ptr(new KDTree(m, this->depth+1));
	    }
	} else {
	    if (this->right != NULL) {
		this->right->add(m);
	    } else {
		this->right = KDTree::Ptr(new KDTree(m, this->depth+1));
	    }
	}
    }
}

int KDTree::count(void) {
    int n = 1;

    if (this->left != NULL)
	n += this->left->count();

    if (this->right != NULL)
	n += this->right->count();

    return n;
}

KDTree::ConstPtr KDTree::findSource(Source const& s) const {

    lsst::afw::geom::Angle ra  = s.getRa();
    lsst::afw::geom::Angle dec = s.getDec();

    lsst::afw::geom::Angle val;
    if (this->axis == 0)
	val = ra;
    else
	val = dec;

    lsst::afw::coord::Coord coord(s.getRa(), s.getDec());
    for (size_t i = 0; i < this->set.size(); i++) {
        // Previous code compared x,y, but those aren't available always now so using RA,Dec.
        // Is this too slow?
        if (coord.angularSeparation(set[i]->getSky()).asArcseconds() < 0.01) {
	    return shared_from_this();
	}
    }

    if (val < this->location[this->axis]) {
	if (this->left != NULL) {
	    return this->left->findSource(s);
	} else {
	    return KDTree::Ptr();
	}
    } else {
	if (this->right != NULL) {
	    return this->right->findSource(s);
	} else {
	    return KDTree::Ptr();
	}
    }

}

KDTree::Ptr KDTree::findNearest(Source const& s) {

    if (this->isLeaf()) {
        return shared_from_this();
    }

    KDTree::Ptr leaf;
    lsst::afw::geom::Angle val;
    if (this->axis == 0) {
	val = s.getRa();
    } else {
	val = s.getDec();
    }

    if (val < this->location[this->axis]) {
	if (this->left != NULL) {
	    if (this->left->isLeaf()) {
		leaf = this->left;
	    } else {
		leaf = this->left->findNearest(s);
	    }
	} else {
	    leaf = this->right->findNearest(s);
	}
	double d_leaf = leaf->distance(s);
	lsst::afw::coord::Coord c;
	if (this->axis == 0) {
	    c = lsst::afw::coord::Coord(val, this->location[1]);
	} else {
	    c = lsst::afw::coord::Coord(this->location[0], val);
	}
	double d_this = this->distance(Source(c));
	if (d_leaf > d_this && this->right != NULL) {
	    KDTree::Ptr leaf2 = this->right->findNearest(s);
	    double d_leaf2 = leaf2->distance(s);
	    if (d_leaf > d_leaf2) {
		leaf = leaf2;
	    }
	}
    } else {
	if (this->right != NULL) {
	    if (this->right->isLeaf()) {
		leaf = this->right;
	    } else {
		leaf = this->right->findNearest(s);
	    }
	} else {
	    leaf = this->left->findNearest(s);
	}
	double d_leaf = leaf->distance(s);
	lsst::afw::coord::Coord c;
	if (this->axis == 0) {
	    c = lsst::afw::coord::Coord(val, this->location[1]);
	} else {
	    c = lsst::afw::coord::Coord(this->location[0], val);
	}
	double d_this = this->distance(Source(c));
	if (d_leaf > d_this && this->left != NULL) {
	    KDTree::Ptr leaf2 = this->left->findNearest(s);
	    double d_leaf2 = leaf2->distance(s);
	    if (d_leaf > d_leaf2) {
		leaf = leaf2;
	    }
	}
    }

    double d_leaf = leaf->distance(s);
    double d_this = this->distance(s);

    if (d_leaf < d_this) {
	return leaf;
    } else {
	return shared_from_this();
    }
}

void KDTree::add(PTR(Source) s, lsst::afw::geom::Angle d_lim) {
    lsst::afw::geom::Angle ra  = s->getRa();
    lsst::afw::geom::Angle dec = s->getDec();

    if (d_lim <= 0) {
        for (size_t i = 0; i < this->set.size(); i++) {
            if (fabs(this->set[i]->getRa()  - ra)  < d_lim &&
                fabs(this->set[i]->getDec() - dec) < d_lim) {
                this->set.push_back(s);
                return;
            }
        }
    }

    double val;
    if (this->axis == 0)
	val = ra;
    else
	val = dec;

    if (val < this->location[this->axis]) {
	if (this->left != NULL) {
	    this->left->add(s, d_lim);
	} else {
	    this->left = KDTree::Ptr(new KDTree(s, this->depth+1));
	}
    } else {
	if (this->right != NULL) {
	    this->right->add(s, d_lim);
	} else {
	    this->right = KDTree::Ptr(new KDTree(s, this->depth+1));
	}
    }
}

SourceGroup KDTree::mergeMat() const {
    SourceGroup sg;
    sg.push_back(this->set);

    if (this->left != NULL) {
        SourceGroup sg_left = this->left->mergeMat();
        for (size_t i = 0; i < sg_left.size(); i++) {
            sg.push_back(sg_left[i]);
        }
    }
    if (this->right != NULL) {
        SourceGroup sg_right = this->right->mergeMat();
        for (size_t i = 0; i < sg_right.size(); i++) {
            sg.push_back(sg_right[i]);
        }
    }

    return sg;
}

SourceGroup KDTree::mergeSource(unsigned int minNumMatch) {
    SourceGroup sg;
    if (this->set.size() >= minNumMatch) {
	double sr = 0.0;
	double sd = 0.0;
	double sm = 0.0;
	double sn = 0.0;
	for (size_t i = 0; i < set.size(); i++) {
	    sr += set[i]->getRa().asDegrees();
	    sd += set[i]->getDec().asDegrees();
	    sm += set[i]->getFlux();
	    sn += 1.0;
	}
	double ra  = sr / sn;
	double dec = sd / sn;
	double mag = sm / sn;
        PTR(Source) source(new Source(lsst::afw::coord::Coord(lsst::afw::geom::Point2D(ra, dec),
                                                              lsst::afw::geom::degrees), mag));
        this->set.insert(set.begin(), source);
	sg.push_back(this->set);
    }

    if (this->left != NULL) {
	SourceGroup sg_left = this->left->mergeSource(minNumMatch);
	for (size_t i = 0; i < sg_left.size(); i++) {
	    sg.push_back(sg_left[i]);
	}
    }
    if (this->right != NULL) {
	SourceGroup sg_right = this->right->mergeSource(minNumMatch);
	for (size_t i = 0; i < sg_right.size(); i++) {
	    sg.push_back(sg_right[i]);
	}
    }

    return sg;
}

void KDTree::printMat() const {
    double ra = set[0]->getRa().asDegrees();
    double dec = set[0]->getDec().asDegrees();

    std::cout << "circle(" << ra << "," << dec << ",5.0\") # color=magenta" << std::endl;

    if (this->left != NULL) {
	this->left->printMat();
    }
    if (this->right != NULL) {
	this->right->printMat();
    }
}

void KDTree::printSource() const {
    double sr = 0.0;
    double sd = 0.0;
    double sn = 0.0;
    for (size_t i = 0; i < set.size(); i++) {
	sr += set[i]->getRa().asDegrees();
	sd += set[i]->getDec().asDegrees();
	sn += 1.0;
    }
    double ra  = sr / sn;
    double dec = sd / sn;

    if (sn >= 2.0)
	std::cout << "circle(" << ra << "," << dec << ",5.0\") # color=red" << std::endl;
    else
	std::cout << "circle(" << ra << "," << dec << ",5.0\")" << std::endl;

    if (this->left != NULL) {
	this->left->printSource();
    }
    if (this->right != NULL) {
	this->right->printSource();
    }
}

KDTree::Ptr
lsst::meas::mosaic::kdtreeMat(SourceMatchGroup &matchList) {

    KDTree::Ptr root = KDTree::Ptr(new KDTree(matchList[0], 0));
    //std::cout << "root->count() : " << root->count() << std::endl;

    for (unsigned int j = 1; j < matchList.size(); j++) {
	for (unsigned int i = 0; i < matchList[j].size(); i++) {
	    root->add(matchList[j][i]);
	}
	//std::cout << "root->count() : " << root->count() << std::endl;
    }

    //std::cout << root->count() << std::endl;

    return root;
}
#if 0
KDTree::Ptr
lsst::meas::mosaic::kdtreeSource(SourceGroup const &sourceSet,
				 KDTree::Ptr rootMat,
				 CcdSet &ccdSet,
				 lsst::afw::geom::Angle d_lim, unsigned int nbrightest) {
    int nchip = ccdSet.size();
    double fluxlim[sourceSet.size()*nchip];

    for (size_t j = 0; j < sourceSet.size(); j++) {
	int k = 0;
	for (CcdSet::iterator it = ccdSet.begin(); it != ccdSet.end(); it++, k++) {
	    std::vector<double> v;
	    for (size_t i = 0; i < sourceSet[j].size(); i++) {
		if (sourceSet[j][i]->getChip() == it->first) {
		    v.push_back(sourceSet[j][i]->getFlux());
		}
	    }
	    if (nbrightest < v.size()) {
		std::sort(v.begin(), v.end(), std::greater<double>());
		fluxlim[j*nchip+k] = v[nbrightest-1];
	    } else {
		fluxlim[j*nchip+k] = 0.0;
	    }
	}
    }
    //std::cout << "(1) " << sourceSet[0].size() << std::endl;

    SourceSet set;
    for (size_t i = 0; i < sourceSet[0].size(); i++) {
	//int k = sourceSet[0][i]->getChip();
	int k = 0;
	for (CcdSet::iterator it = ccdSet.begin(); it != ccdSet.end(); it++, k++) {
	    if (sourceSet[0][i]->getChip() == it->first) break;
	}
        if (sourceSet[0][i]->getFlux() >= fluxlim[k] &&
	    rootMat->findSource(*sourceSet[0][i]) == NULL) {
	    set.push_back(sourceSet[0][i]);
	}
    }
    //std::cout << "(2) " << set.size() << std::endl;

    KDTree::Ptr rootSource;
    if (set.size() > 0) {
        rootSource = KDTree::Ptr(new KDTree(set, 0));
    }

    //std::cout << "(3) " << rootSource->count() << std::endl;
    for (size_t j = 1; j < sourceSet.size(); j++) {
	for (size_t i = 0; i < sourceSet[j].size(); i++) {
	    //int k = sourceSet[j][i]->getChip();
	    int k = 0;
	    for (CcdSet::iterator it = ccdSet.begin(); it != ccdSet.end(); it++, k++) {
		if (sourceSet[j][i]->getChip() == it->first) break;
	    }
	    if (sourceSet[j][i]->getFlux() >= fluxlim[j*nchip+k] &&
		rootMat->findSource(*sourceSet[j][i]) == NULL) {
                if (rootSource) {
                    KDTree::Ptr leaf = rootSource->findNearest(*sourceSet[j][i]);
                    if (leaf->distance(*sourceSet[j][i]) < d_lim) {
                        leaf->set.push_back(sourceSet[j][i]);
                    } else {
                        rootSource->add(sourceSet[j][i]);
                    }
                } else {
                    rootSource = KDTree::Ptr(new KDTree(sourceSet[j][i], 0));
                }
	    }
	}
	//std::cout << "(3) " << rootSource->count() << std::endl;
    }
    //std::cout << "(4) " << rootSource->count() << std::endl;

    return rootSource;
}
#endif
KDTree::Ptr
lsst::meas::mosaic::kdtreeSource(SourceGroup const &sourceSet,
				 KDTree::Ptr rootMat,
				 CcdSet &ccdSet,
				 lsst::afw::geom::Angle d_lim) {
    //std::cout << "(1) " << sourceSet[0].size() << std::endl;

    SourceSet set;
    for (size_t i = 0; i < sourceSet[0].size(); i++) {
        if (rootMat->findSource(*sourceSet[0][i]) == NULL) {
	    set.push_back(sourceSet[0][i]);
	}
    }
    //std::cout << "(2) " << set.size() << std::endl;

    KDTree::Ptr rootSource;
    if (set.size() > 0) {
        rootSource = KDTree::Ptr(new KDTree(set, 0));
    }

    //std::cout << "(3) " << rootSource->count() << std::endl;
    for (size_t j = 1; j < sourceSet.size(); j++) {
	for (size_t i = 0; i < sourceSet[j].size(); i++) {
	    if (rootMat->findSource(*sourceSet[j][i]) == NULL) {
                if (rootSource) {
                    KDTree::Ptr leaf = rootSource->findNearest(*sourceSet[j][i]);
                    if (leaf->distance(*sourceSet[j][i]) < d_lim.asDegrees()) {
                        leaf->set.push_back(sourceSet[j][i]);
                    } else {
                        rootSource->add(sourceSet[j][i]);
                    }
                } else {
                    rootSource = KDTree::Ptr(new KDTree(sourceSet[j][i], 0));
                }
	    }
	}
	//std::cout << "(3) " << rootSource->count() << std::endl;
    }
    //std::cout << "(4) " << rootSource->count() << std::endl;

    return rootSource;
}

double calXi(double a, double d, double A, double D) {
    return cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_a(double a, double d, double A, double D) {
    return cos(D)*pow(cos(d),2.)*pow(sin(a-A),2.)/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	  +cos(d)*cos(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_d(double a, double d, double A, double D) {
    return -cos(d)*sin(a-A)*(sin(D)*cos(d)-cos(D)*sin(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   -sin(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_A(double a, double d, double A, double D) {
    return -cos(D)*cos(d)*cos(d)*sin(a-A)*sin(a-A)/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   -cos(d)*cos(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_D(double a, double d, double A, double D) {
    return -cos(d)*sin(a-A)*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.);
}

double calEta(double a, double d, double A, double D) {
    return (cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_a(double a, double d, double A, double D) {
    return cos(D)*cos(d)*sin(a-A)*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	  +sin(D)*cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_d(double a, double d, double A, double D) {
    return -(sin(D)*cos(d)-cos(D)*sin(d)*cos(a-A))*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   +(cos(D)*cos(d)+sin(D)*sin(d)*cos(a-A))/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_A(double a, double d, double A, double D) {
    return -cos(D)*cos(d)*sin(a-A)*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   -sin(D)*cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_D(double a, double d, double A, double D) {
    return -pow(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A),2.)/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)-1.;
}

Eigen::VectorXd solveMatrix_MKL(long size, Eigen::MatrixXd &a_data, Eigen::VectorXd &b_data) {
    //char L = 'L';
    lapack::MKL_INT n = size;
    lapack::MKL_INT nrhs = 1;
    lapack::MKL_INT lda = size;
    std::vector<lapack::MKL_INT> ipiv(size);
    lapack::MKL_INT ldb = size;
    lapack::MKL_INT info = 0;

    lapack::dgesv(&n, &nrhs, &a_data(0), &lda, ipiv.data(), &b_data(0), &ldb, &info);
    //lapack::dposv(&L, &n, &nrhs, &a_data(0), &lda, &b_data(0), &ldb, &info);
    if(info != 0){
	throw std::runtime_error((boost::format("solving linear equation failed: dgesv returned %1%") % info).str());
    }

    Eigen::VectorXd c_data(size);
    for (int i = 0; i < size; i++) {
	c_data(i) = b_data(i);
    }

    return c_data;
}

Eigen::VectorXd solveMatrix_Eigen(long size, Eigen::MatrixXd &a, Eigen::VectorXd &b) {
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(a);
    return lu.solve(b);
}

Eigen::VectorXd solveMatrix(long size, Eigen::MatrixXd &a_data, Eigen::VectorXd &b_data) {
    if(lapack::isLapackAvailable) {
	return solveMatrix_MKL(size, a_data, b_data);
    }
    else {
	return solveMatrix_Eigen(size, a_data, b_data);
    }
}

Eigen::VectorXd solveForCoeff(std::vector<Obs::Ptr>& objList, Poly::Ptr p) {
    int ncoeff = p->ncoeff;
    int size = 2 * ncoeff + 2;

    int *xorder = p->xorder;
    int *yorder = p->yorder;

    Eigen::MatrixXd a_data = Eigen::MatrixXd::Zero(size, size);
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(size);

    Eigen::VectorXd pu(ncoeff);
    Eigen::VectorXd pv(ncoeff);

    for (size_t k = 0; k < objList.size(); k++) {
	Obs::Ptr o = objList[k];
	if (o->good) {
	    for (int j = 0; j < ncoeff; j++) {
		pu(j) = pow(o->u, xorder[j]);
		pv(j) = pow(o->v, yorder[j]);
	    }
	    for (int j = 0; j < ncoeff; j++) {
		b_data(j)        += o->xi  * pu(j) * pv(j);
		b_data(j+ncoeff) += o->eta * pu(j) * pv(j);
		for (int i = 0; i < ncoeff; i++) {
		    a_data(i,        j       ) += pu(j) * pv(j) * pu(i) * pv(i);
		    a_data(i+ncoeff, j+ncoeff) += pu(j) * pv(j) * pu(i) * pv(i);
		}
		a_data(j,          2*ncoeff  ) -= pu(j) * pv(j) * o->xi_A;
		a_data(j,          2*ncoeff+1) -= pu(j) * pv(j) * o->xi_D;
		a_data(j+ncoeff,   2*ncoeff  ) -= pu(j) * pv(j) * o->eta_A;
		a_data(j+ncoeff,   2*ncoeff+1) -= pu(j) * pv(j) * o->eta_D;
		a_data(2*ncoeff,   j         ) -= pu(j) * pv(j) * o->xi_A;
		a_data(2*ncoeff+1, j         ) -= pu(j) * pv(j) * o->xi_D;
		a_data(2*ncoeff,   j+ncoeff  ) -= pu(j) * pv(j) * o->eta_A;
		a_data(2*ncoeff+1, j+ncoeff  ) -= pu(j) * pv(j) * o->eta_D;
	    }
	    a_data(2*ncoeff  , 2*ncoeff  ) += o->xi_A * o->xi_A + o->eta_A * o->eta_A;
	    a_data(2*ncoeff  , 2*ncoeff+1) += o->xi_A * o->xi_D + o->eta_A * o->eta_D;
	    a_data(2*ncoeff+1, 2*ncoeff  ) += o->xi_A * o->xi_D + o->eta_A * o->eta_D;
	    a_data(2*ncoeff+1, 2*ncoeff+1) += o->xi_D * o->xi_D + o->eta_D * o->eta_D;
	    b_data(2*ncoeff)   -= o->xi * o->xi_A + o->eta * o->eta_A;
	    b_data(2*ncoeff+1) -= o->xi * o->xi_D + o->eta * o->eta_D;
	}
    }

    Eigen::VectorXd coeff = solveMatrix(size, a_data, b_data);
    //Eigen::VectorXd coeff = a_data.partialPivLu().solve(b_data);

    return coeff;
}

Eigen::VectorXd solveForCoeffWithOffset(std::vector<Obs::Ptr>& objList, Coeff::Ptr& c, Poly::Ptr p) {
    int ncoeff = p->ncoeff;
    int size = 2 * ncoeff + 2;

    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double *a = c->a;
    double *b = c->b;

    Eigen::MatrixXd a_data = Eigen::MatrixXd::Zero(size, size);
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(size);

    Eigen::VectorXd pu(ncoeff);
    Eigen::VectorXd pv(ncoeff);

    for (size_t i = 0; i < objList.size(); i++) {
	Obs::Ptr o = objList[i];
	if (o->good) {
	    double Ax = o->xi;
	    double Ay = o->eta;
	    double Bx = 0.0;
	    double By = 0.0;
	    double Cx = 0.0;
	    double Cy = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		pu(k) = pow(o->u, xorder[k]);
		pv(k) = pow(o->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[k] * pu(k) * pv(k);
		Ay -= b[k] * pu(k) * pv(k);
		Bx += a[k] * pow(o->u, xorder[k]-1) * pv(k) * xorder[k];
		By += b[k] * pow(o->u, xorder[k]-1) * pv(k) * xorder[k];
		Cx += a[k] * pu(k) * pow(o->v, yorder[k]-1) * yorder[k];
		Cy += b[k] * pu(k) * pow(o->v, yorder[k]-1) * yorder[k];
	    }
	    for (int k = 0; k < ncoeff; k++) {
		b_data(k)        += Ax * pu(k) * pv(k);
		b_data(k+ncoeff) += Ay * pu(k) * pv(k);
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data(j,        k       ) += pu(j) * pv(j) * pu(k) * pv(k);
		    a_data(j+ncoeff, k+ncoeff) += pu(j) * pv(j) * pu(k) * pv(k);
		}

		// coeff x offset
		a_data(k         , ncoeff*2  ) += Bx * pu(k) * pv(k);
		a_data(k         , ncoeff*2+1) += Cx * pu(k) * pv(k);
		a_data(k+ncoeff  , ncoeff*2  ) += By * pu(k) * pv(k);
		a_data(k+ncoeff  , ncoeff*2+1) += Cy * pu(k) * pv(k);
		a_data(ncoeff*2  , k         ) += Bx * pu(k) * pv(k);
		a_data(ncoeff*2+1, k         ) += Cx * pu(k) * pv(k);
		a_data(ncoeff*2  , k+ncoeff  ) += By * pu(k) * pv(k);
		a_data(ncoeff*2+1, k+ncoeff  ) += Cy * pu(k) * pv(k);
	    }

	    // offset x offset
	    a_data(ncoeff*2  , ncoeff*2  ) += Bx * Bx + By * By;
	    a_data(ncoeff*2  , ncoeff*2+1) += Bx * Cx + By * Cy;
	    a_data(ncoeff*2+1, ncoeff*2  ) += Cx * Bx + Cy * By;
	    a_data(ncoeff*2+1, ncoeff*2+1) += Cx * Cx + Cy * Cy;

	    b_data(ncoeff*2  ) += Ax * Bx + Ay * By;
	    b_data(ncoeff*2+1) += Ax * Cx + Ay * Cy;
	}
    }

    Eigen::VectorXd coeff = solveMatrix(size, a_data, b_data);
//    Eigen::VectorXd coeff = a_data.partialPivLu().solve(b_data);

    return coeff;
}

double calcChi(std::vector<Obs::Ptr>& objList, Eigen::VectorXd &a, Poly::Ptr p) {
    int ncoeff = p->ncoeff;

    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double chi2 = 0.0;
    for (size_t k = 0; k < objList.size(); k++) {
	Obs::Ptr o = objList[k];
	if (o->good) {
	    double Ax = o->xi;
	    double Ay = o->eta;
	    for (int i = 0; i < ncoeff; i++) {
		Ax -= a(i)        * pow(o->u, xorder[i]) * pow(o->v, yorder[i]);
		Ay -= a(i+ncoeff) * pow(o->u, xorder[i]) * pow(o->v, yorder[i]);
	    }
	    Ax += (o->xi_A  * a(2*ncoeff) + o->xi_D  * a(2*ncoeff+1));
	    Ay += (o->eta_A * a(2*ncoeff) + o->eta_D * a(2*ncoeff+1));
	    chi2 += Ax * Ax + Ay * Ay;
	}
    }

    return chi2;
}

double flagObj(std::vector<Obs::Ptr>& objList, Eigen::VectorXd &a, Poly::Ptr p, double e2) {
    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double chi2 = 0.0;
    int nrejected = 0;
    for (size_t j = 0; j < objList.size(); j++) {
	Obs::Ptr o = objList[j];
	double Ax = 0.0;
	double Ay = 0.0;
	for (int i = 0; i < ncoeff; i++) {
	    double pu = pow(o->u, xorder[i]);
	    double pv = pow(o->v, yorder[i]);
	    Ax += a(i)        * pu * pv;
	    Ay += a(i+ncoeff) * pu * pv;
	}
	Ax -= (o->xi_A  * a(2*ncoeff) + o->xi_D  * a(2*ncoeff+1));
	Ay -= (o->eta_A * a(2*ncoeff) + o->eta_D * a(2*ncoeff+1));
	double r2 = pow(o->xi - Ax, 2) + pow(o->eta - Ay, 2);
	if (r2 > e2) {
	    o->good = false;
	    nrejected++;
	} else {
	    //o->good = true;
	}
    }
    printf("nrejected = %d\n", nrejected);

    return chi2;
}

Eigen::VectorXd
solveLinApprox(std::vector<Obs::Ptr>& o, CoeffSet& coeffVec, int nchip, Poly::Ptr p,
	       bool solveCcd=true,
	       bool allowRotation=true,
	       double catRMS=0.0)
{
    int nobs  = o.size();
    int nexp = coeffVec.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    std::map<int, double*> a;
    std::map<int, double*> b;
    for (CoeffSet::iterator it = coeffVec.begin(); it != coeffVec.end(); it++) {
	a.insert(std::map<int, double*>::value_type(it->first, it->second->a));
	b.insert(std::map<int, double*>::value_type(it->first, it->second->b));
    }

    long size, np = 0;
    if (solveCcd) {
	if (allowRotation) {
	    size = 2 * ncoeff * nexp + 3 * nchip + 1;
	    np = 3;
	} else {
	    size = 2 * ncoeff * nexp + 2 * nchip;
	    np = 2;
	}
    } else {
	size = 2 * ncoeff * nexp;
    }
    std::cout << "size: " << size << std::endl;

    Eigen::MatrixXd a_data = Eigen::MatrixXd::Zero(size, size);
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(size);

    Eigen::VectorXd pu(ncoeff);
    Eigen::VectorXd pv(ncoeff);

    double isx2 = 1.0;
    double isy2 = 1.0;

    if (solveCcd) {
	for (int i = 0; i < nobs; i++) {
	    if (!o[i]->good) continue;
	    double Ax = o[i]->xi;
	    double Ay = o[i]->eta;
	    double Bx = 0.0;
	    double By = 0.0;
	    double Cx = 0.0;
	    double Cy = 0.0;
	    double Dx = 0.0;
	    double Dy = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		pu(k) = pow(o[i]->u, xorder[k]);
		pv(k) = pow(o[i]->v, yorder[k]);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu(k)   * pv(k);
		Ay -= b[o[i]->iexp][k] * pu(k)   * pv(k);
		Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		Cx += a[o[i]->iexp][k] * pu(k)   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[o[i]->iexp][k] * pu(k)   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Dx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
		Dy += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
	    }
	    double dxi  = Bx * o[i]->xerr + Cx * o[i]->yerr;
	    double deta = By * o[i]->xerr + Cy * o[i]->yerr;
	    isx2 = 1.0 / (pow(dxi,  2) + pow(catRMS, 2));
	    isy2 = 1.0 / (pow(deta, 2) + pow(catRMS, 2));

	    for (int k = 0; k < ncoeff; k++) {
		b_data(k+       ncoeff*2*o[i]->jexp) += Ax * pu(k) * pv(k) * isx2;
		b_data(k+ncoeff+ncoeff*2*o[i]->jexp) += Ay * pu(k) * pv(k) * isy2;
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data(j+       ncoeff*2*o[i]->jexp, k+       ncoeff*2*o[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isx2;
		    a_data(j+ncoeff+ncoeff*2*o[i]->jexp, k+ncoeff+ncoeff*2*o[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isy2;
		}

		// coeff x chip
		a_data(k+       ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np  ) += Bx * pu(k) * pv(k) * isx2;
		a_data(k+       ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np+1) += Cx * pu(k) * pv(k) * isx2;
		a_data(k+ncoeff+ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np  ) += By * pu(k) * pv(k) * isy2;
		a_data(k+ncoeff+ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np+1) += Cy * pu(k) * pv(k) * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np  , k+       ncoeff*2*o[i]->jexp) += Bx * pu(k) * pv(k) * isx2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+1, k+       ncoeff*2*o[i]->jexp) += Cx * pu(k) * pv(k) * isx2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np  , k+ncoeff+ncoeff*2*o[i]->jexp) += By * pu(k) * pv(k) * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+1, k+ncoeff+ncoeff*2*o[i]->jexp) += Cy * pu(k) * pv(k) * isy2;
		if (allowRotation) {
		    a_data(k+       ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np+2) += Dx * pu(k) * pv(k) * isx2;
		    a_data(k+ncoeff+ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np+2) += Dy * pu(k) * pv(k) * isy2;
		    a_data(ncoeff*2*nexp+o[i]->jchip*np+2, k+       ncoeff*2*o[i]->jexp) += Dx * pu(k) * pv(k) * isx2;
		    a_data(ncoeff*2*nexp+o[i]->jchip*np+2, k+ncoeff+ncoeff*2*o[i]->jexp) += Dy * pu(k) * pv(k) * isy2;
		}
	    }

	    // chip x chip
	    a_data(ncoeff*2*nexp+o[i]->jchip*np  , ncoeff*2*nexp+o[i]->jchip*np  ) += Bx * Bx * isx2 + By * By * isy2;
	    a_data(ncoeff*2*nexp+o[i]->jchip*np  , ncoeff*2*nexp+o[i]->jchip*np+1) += Bx * Cx * isx2 + By * Cy * isy2;
	    a_data(ncoeff*2*nexp+o[i]->jchip*np+1, ncoeff*2*nexp+o[i]->jchip*np  ) += Cx * Bx * isx2 + Cy * By * isy2;
	    a_data(ncoeff*2*nexp+o[i]->jchip*np+1, ncoeff*2*nexp+o[i]->jchip*np+1) += Cx * Cx * isx2 + Cy * Cy * isy2;
	    if (allowRotation) {
		a_data(ncoeff*2*nexp+o[i]->jchip*np  , ncoeff*2*nexp+o[i]->jchip*np+2) += Bx * Dx * isx2 + By * Dy * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+1, ncoeff*2*nexp+o[i]->jchip*np+2) += Cx * Dx * isx2 + Cy * Dy * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+2, ncoeff*2*nexp+o[i]->jchip*np  ) += Dx * Bx * isx2 + Dy * By * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+2, ncoeff*2*nexp+o[i]->jchip*np+1) += Dx * Cx * isx2 + Dy * Cy * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+2, ncoeff*2*nexp+o[i]->jchip*np+2) += Dx * Dx * isx2 + Dy * Dy * isy2;
	    }

	    b_data[ncoeff*2*nexp+o[i]->jchip*np  ] += Ax * Bx * isx2 + Ay * By * isy2;
	    b_data[ncoeff*2*nexp+o[i]->jchip*np+1] += Ax * Cx * isx2 + Ay * Cy * isy2;
	    if (allowRotation) {
		b_data[ncoeff*2*nexp+o[i]->jchip*np+2] += Ax * Dx * isx2 + Ay * Dy * isy2;
	    }
	}

	if (allowRotation) {
	    // \Sum d_theta = 0.0
	    for (int i = 0; i < nchip; i++) {
		a_data(ncoeff*2*nexp+i*np+2, ncoeff*2*nexp+nchip*np) = 1;
		a_data(ncoeff*2*nexp+nchip*np, ncoeff*2*nexp+i*np+2) = 1;
	    }
	}
    } else {
	for (int i = 0; i < nobs; i++) {
	    if (!o[i]->good) continue;
	    double Ax = o[i]->xi;
	    double Ay = o[i]->eta;
	    double Bx = 0.0;
	    double By = 0.0;
	    double Cx = 0.0;
	    double Cy = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		pu(k) = pow(o[i]->u, xorder[k]);
		pv(k) = pow(o[i]->v, yorder[k]);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu(k)   * pv(k);
		Ay -= b[o[i]->iexp][k] * pu(k)   * pv(k);
		Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		Cx += a[o[i]->iexp][k] * pu(k)   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[o[i]->iexp][k] * pu(k)   * pow(o[i]->v, yorder[k]-1) * yorder[k];
	    }
	    double dxi  = Bx * o[i]->xerr + Cx * o[i]->yerr;
	    double deta = By * o[i]->xerr + Cy * o[i]->yerr;
	    isx2 = 1.0 / (pow(dxi,  2) + pow(catRMS, 2));
	    isy2 = 1.0 / (pow(deta, 2) + pow(catRMS, 2));

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*o[i]->jexp] += Ax * pu(k) * pv(k) * isx2;
		b_data[k+ncoeff+ncoeff*2*o[i]->jexp] += Ay * pu(k) * pv(k) * isy2;
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data(j+       ncoeff*2*o[i]->jexp, k+       ncoeff*2*o[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isx2;
		    a_data(j+ncoeff+ncoeff*2*o[i]->jexp, k+ncoeff+ncoeff*2*o[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isy2;
		}
	    }
	}
    }

    Eigen::VectorXd coeff = solveMatrix(size, a_data, b_data);
//    Eigen::VectorXd coeff = a_data.partialPivLu().solve(b_data);

    return coeff;
}

Eigen::VectorXd
solveLinApprox_Star(std::vector<Obs::Ptr>& o, std::vector<Obs::Ptr>& s, int nstar,
		    CoeffSet coeffVec, int nchip, Poly::Ptr p,
		    bool solveCcd=true,
		    bool allowRotation=true,
		    double catRMS=0.0)
{
    int nobs  = o.size();
    int nSobs = s.size();
    int nexp = coeffVec.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    std::map<int, double*> a;
    std::map<int, double*> b;
    for (CoeffSet::iterator it = coeffVec.begin(); it != coeffVec.end(); it++) {
	a.insert(std::map<int, double*>::value_type(it->first, it->second->a));
	b.insert(std::map<int, double*>::value_type(it->first, it->second->b));
    }

    int* num = new int[nstar];
    for (int i = 0; i < nstar; i++) {
	num[i] = 0;
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good) {
	    num[s[i]->istar] += 1;
	}
    }
    std::vector<int> v_istar;
    for (int i = 0; i < nstar; i++) {
	if (num[i] >= 2) {
	    v_istar.push_back(i);
	}
    }
    delete [] num;
    int nstar2 = v_istar.size();
    std::cout << "nstar: " << nstar2 << std::endl;

    for (int i = 0; i < nSobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), s[i]->istar);
	if (it != v_istar.end()) {
	    s[i]->jstar = it - v_istar.begin();
	} else {
	    s[i]->jstar = -1;
	}
    }

    long size, size0, np = 0;
    if (solveCcd) {
	if (allowRotation) {
	    size  = 2 * ncoeff * nexp + 3 * nchip + 1 + nstar2 * 2;
	    size0 = 2 * ncoeff * nexp + 3 * nchip + 1;
	    np = 3;
	} else {
	    size  = 2 * ncoeff * nexp + 2 * nchip + nstar2 * 2;
	    size0 = 2 * ncoeff * nexp + 2 * nchip;
	    np = 2;
	}
    } else {
	size  = 2 * ncoeff * nexp + nstar2 * 2;
	size0 = 2 * ncoeff * nexp;
    }

    std::cout << "size : " << size << std::endl;

    Eigen::MatrixXd a_data;
    try {
	a_data = Eigen::MatrixXd::Zero(size, size);
    } catch (std::bad_alloc) {
	std::cerr << "Memory allocation error: for a_data" << std::endl;
	fprintf(stderr, "You need %5.1f GB memory\n", size*size*sizeof(double)/double(1024*1024*1024));
	abort();
    }
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(size);

    Eigen::VectorXd pu(ncoeff);
    Eigen::VectorXd pv(ncoeff);

    int numObsGood = 0, numStarGood = 0;

    double isx2 = 1.0;
    double isy2 = 1.0;

    if (solveCcd) {
	for (int i = 0; i < nobs; i++) {
	    if (!o[i]->good) continue;
            ++numObsGood;
	    double Ax = o[i]->xi;
	    double Ay = o[i]->eta;
	    double Bx = 0.0;
	    double By = 0.0;
	    double Cx = 0.0;
	    double Cy = 0.0;
	    double Dx = 0.0;
	    double Dy = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		pu(k) = pow(o[i]->u, xorder[k]);
		pv(k) = pow(o[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu(k)   * pv(k);
		Ay -= b[o[i]->iexp][k] * pu(k)   * pv(k);
		Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		Cx += a[o[i]->iexp][k] * pu(k)   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[o[i]->iexp][k] * pu(k)   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Dx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
		Dy += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
	    }
	    double dxi  = Bx * o[i]->xerr + Cx * o[i]->yerr;
	    double deta = By * o[i]->xerr + Cy * o[i]->yerr;
	    isx2 = 1.0 / (pow(dxi,  2) + pow(catRMS, 2));
	    isy2 = 1.0 / (pow(deta, 2) + pow(catRMS, 2));

	    for (int k = 0; k < ncoeff; k++) {
		b_data(k+       ncoeff*2*o[i]->jexp) += Ax * pu(k) * pv(k) * isx2;
		b_data(k+ncoeff+ncoeff*2*o[i]->jexp) += Ay * pu(k) * pv(k) * isy2;
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data(j+       ncoeff*2*o[i]->jexp, k+       ncoeff*2*o[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isx2;
		    a_data(j+ncoeff+ncoeff*2*o[i]->jexp, k+ncoeff+ncoeff*2*o[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isy2;
		}

		// coeff x chip
		a_data(k+       ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np  ) += Bx * pu(k) * pv(k) * isx2;
		a_data(k+       ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np+1) += Cx * pu(k) * pv(k) * isx2;
		a_data(k+ncoeff+ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np  ) += By * pu(k) * pv(k) * isy2;
		a_data(k+ncoeff+ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np+1) += Cy * pu(k) * pv(k) * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np  , k+       ncoeff*2*o[i]->jexp) += Bx * pu(k) * pv(k) * isx2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+1, k+       ncoeff*2*o[i]->jexp) += Cx * pu(k) * pv(k) * isx2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np  , k+ncoeff+ncoeff*2*o[i]->jexp) += By * pu(k) * pv(k) * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+1, k+ncoeff+ncoeff*2*o[i]->jexp) += Cy * pu(k) * pv(k) * isy2;
		if (allowRotation) {
		    a_data(k+       ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np+2) += Dx * pu(k) * pv(k) * isx2;
		    a_data(k+ncoeff+ncoeff*2*o[i]->jexp, ncoeff*2*nexp+o[i]->jchip*np+2) += Dy * pu(k) * pv(k) * isy2;
		    a_data(ncoeff*2*nexp+o[i]->jchip*np+2, k+       ncoeff*2*o[i]->jexp) += Dx * pu(k) * pv(k) * isx2;
		    a_data(ncoeff*2*nexp+o[i]->jchip*np+2, k+ncoeff+ncoeff*2*o[i]->jexp) += Dy * pu(k) * pv(k) * isy2;
		}
	    }

	    // chip x chip
	    a_data(ncoeff*2*nexp+o[i]->jchip*np  , ncoeff*2*nexp+o[i]->jchip*np  ) += Bx * Bx * isx2 + By * By * isy2;
	    a_data(ncoeff*2*nexp+o[i]->jchip*np  , ncoeff*2*nexp+o[i]->jchip*np+1) += Bx * Cx * isx2 + By * Cy * isy2;
	    a_data(ncoeff*2*nexp+o[i]->jchip*np+1, ncoeff*2*nexp+o[i]->jchip*np  ) += Cx * Bx * isx2 + Cy * By * isy2;
	    a_data(ncoeff*2*nexp+o[i]->jchip*np+1, ncoeff*2*nexp+o[i]->jchip*np+1) += Cx * Cx * isx2 + Cy * Cy * isy2;
	    if (allowRotation) {
		a_data(ncoeff*2*nexp+o[i]->jchip*np  , ncoeff*2*nexp+o[i]->jchip*np+2) += Bx * Dx * isx2 + By * Dy * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+1, ncoeff*2*nexp+o[i]->jchip*np+2) += Cx * Dx * isx2 + Cy * Dy * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+2, ncoeff*2*nexp+o[i]->jchip*np  ) += Dx * Bx * isx2 + Dy * By * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+2, ncoeff*2*nexp+o[i]->jchip*np+1) += Dx * Cx * isx2 + Dy * Cy * isy2;
		a_data(ncoeff*2*nexp+o[i]->jchip*np+2, ncoeff*2*nexp+o[i]->jchip*np+2) += Dx * Dx * isx2 + Dy * Dy * isy2;
	    }

	    b_data(ncoeff*2*nexp+o[i]->jchip*np  ) += Ax * Bx * isx2 + Ay * By * isy2;
	    b_data(ncoeff*2*nexp+o[i]->jchip*np+1) += Ax * Cx * isx2 + Ay * Cy * isy2;
	    if (allowRotation) {
		b_data(ncoeff*2*nexp+o[i]->jchip*np+2) += Ax * Dx * isx2 + Ay * Dy * isy2;
	    }
	}

	for (int i = 0; i < nSobs; i++) {
	    if (!s[i]->good || s[i]->jstar == -1) continue;
            ++numStarGood;
	    double Ax = s[i]->xi;
	    double Ay = s[i]->eta;
	    double Bx = 0.0;
	    double By = 0.0;
	    double Cx = 0.0;
	    double Cy = 0.0;
	    double Dx = 0.0;
	    double Dy = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		pu(k) = pow(s[i]->u, xorder[k]);
		pv(k) = pow(s[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[s[i]->iexp][k] * pu(k)   * pv(k);
		Ay -= b[s[i]->iexp][k] * pu(k)   * pv(k);
		Bx += a[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		By += b[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		Cx += a[s[i]->iexp][k] * pu(k)   * pow(s[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[s[i]->iexp][k] * pu(k)   * pow(s[i]->v, yorder[k]-1) * yorder[k];
		Dx += a[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pow(s[i]->v, yorder[k]-1) * (-xorder[k]*s[i]->v*s[i]->v0+yorder[k]*s[i]->u*s[i]->u0);
		Dy += b[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pow(s[i]->v, yorder[k]-1) * (-xorder[k]*s[i]->v*s[i]->v0+yorder[k]*s[i]->u*s[i]->u0);
	    }
	    double dxi  = Bx * s[i]->xerr + Cx * s[i]->yerr;
	    double deta = By * s[i]->xerr + Cy * s[i]->yerr;
	    isx2 = 1.0 / pow(dxi,  2);
	    isy2 = 1.0 / pow(deta, 2);

	    for (int k = 0; k < ncoeff; k++) {
		b_data(k+       ncoeff*2*s[i]->jexp) += Ax * pu(k) * pv(k) * isx2;
		b_data(k+ncoeff+ncoeff*2*s[i]->jexp) += Ay * pu(k) * pv(k) * isy2;
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data(j+       ncoeff*2*s[i]->jexp, k+       ncoeff*2*s[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isx2;
		    a_data(j+ncoeff+ncoeff*2*s[i]->jexp, k+ncoeff+ncoeff*2*s[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isy2;
		}

		// coeff x chip
		a_data(k+       ncoeff*2*s[i]->jexp, ncoeff*2*nexp+s[i]->jchip*np  ) += Bx * pu(k) * pv(k) * isx2;
		a_data(k+       ncoeff*2*s[i]->jexp, ncoeff*2*nexp+s[i]->jchip*np+1) += Cx * pu(k) * pv(k) * isx2;
		a_data(k+ncoeff+ncoeff*2*s[i]->jexp, ncoeff*2*nexp+s[i]->jchip*np  ) += By * pu(k) * pv(k) * isy2;
		a_data(k+ncoeff+ncoeff*2*s[i]->jexp, ncoeff*2*nexp+s[i]->jchip*np+1) += Cy * pu(k) * pv(k) * isy2;
		a_data(ncoeff*2*nexp+s[i]->jchip*np  , k+       ncoeff*2*s[i]->jexp) += Bx * pu(k) * pv(k) * isx2;
		a_data(ncoeff*2*nexp+s[i]->jchip*np+1, k+       ncoeff*2*s[i]->jexp) += Cx * pu(k) * pv(k) * isx2;
		a_data(ncoeff*2*nexp+s[i]->jchip*np  , k+ncoeff+ncoeff*2*s[i]->jexp) += By * pu(k) * pv(k) * isy2;
		a_data(ncoeff*2*nexp+s[i]->jchip*np+1, k+ncoeff+ncoeff*2*s[i]->jexp) += Cy * pu(k) * pv(k) * isy2;
		if (allowRotation) {
		    a_data(k+       ncoeff*2*s[i]->jexp, ncoeff*2*nexp+s[i]->jchip*np+2) += Dx * pu(k) * pv(k) * isx2;
		    a_data(k+ncoeff+ncoeff*2*s[i]->jexp, ncoeff*2*nexp+s[i]->jchip*np+2) += Dy * pu(k) * pv(k) * isy2;
		    a_data(ncoeff*2*nexp+s[i]->jchip*np+2, k+       ncoeff*2*s[i]->jexp) += Dx * pu(k) * pv(k) * isx2;
		    a_data(ncoeff*2*nexp+s[i]->jchip*np+2, k+ncoeff+ncoeff*2*s[i]->jexp) += Dy * pu(k) * pv(k) * isy2;
		}

		// coeff x star
		a_data(k+       ncoeff*2*s[i]->jexp, size0+s[i]->jstar*2  ) -= s[i]->xi_a  * pu(k) * pv(k) * isx2;
		a_data(k+       ncoeff*2*s[i]->jexp, size0+s[i]->jstar*2+1) -= s[i]->xi_d  * pu(k) * pv(k) * isx2;
		a_data(k+ncoeff+ncoeff*2*s[i]->jexp, size0+s[i]->jstar*2  ) -= s[i]->eta_a * pu(k) * pv(k) * isy2;
		a_data(k+ncoeff+ncoeff*2*s[i]->jexp, size0+s[i]->jstar*2+1) -= s[i]->eta_d * pu(k) * pv(k) * isy2;
		a_data(size0+s[i]->jstar*2  , k+       ncoeff*2*s[i]->jexp) -= s[i]->xi_a  * pu(k) * pv(k) * isx2;
		a_data(size0+s[i]->jstar*2+1, k+       ncoeff*2*s[i]->jexp) -= s[i]->xi_d  * pu(k) * pv(k) * isx2;
		a_data(size0+s[i]->jstar*2  , k+ncoeff+ncoeff*2*s[i]->jexp) -= s[i]->eta_a * pu(k) * pv(k) * isy2;
		a_data(size0+s[i]->jstar*2+1, k+ncoeff+ncoeff*2*s[i]->jexp) -= s[i]->eta_d * pu(k) * pv(k) * isy2;
	    }

	    // chip x chip
	    a_data(ncoeff*2*nexp+s[i]->jchip*np  , ncoeff*2*nexp+s[i]->jchip*np  ) += Bx * Bx * isx2 + By * By * isy2;
	    a_data(ncoeff*2*nexp+s[i]->jchip*np  , ncoeff*2*nexp+s[i]->jchip*np+1) += Bx * Cx * isx2 + By * Cy * isy2;
	    a_data(ncoeff*2*nexp+s[i]->jchip*np+1, ncoeff*2*nexp+s[i]->jchip*np  ) += Cx * Bx * isx2 + Cy * By * isy2;
	    a_data(ncoeff*2*nexp+s[i]->jchip*np+1, ncoeff*2*nexp+s[i]->jchip*np+1) += Cx * Cx * isx2 + Cy * Cy * isy2;
	    if (allowRotation) {
		a_data(ncoeff*2*nexp+s[i]->jchip*np  , ncoeff*2*nexp+s[i]->jchip*np+2) += Bx * Dx * isx2 + By * Dy * isy2;
		a_data(ncoeff*2*nexp+s[i]->jchip*np+1, ncoeff*2*nexp+s[i]->jchip*np+2) += Cx * Dx * isx2 + Cy * Dy * isy2;
		a_data(ncoeff*2*nexp+s[i]->jchip*np+2, ncoeff*2*nexp+s[i]->jchip*np  ) += Dx * Bx * isx2 + Dy * By * isy2;
		a_data(ncoeff*2*nexp+s[i]->jchip*np+2, ncoeff*2*nexp+s[i]->jchip*np+1) += Dx * Cx * isx2 + Dy * Cy * isy2;
		a_data(ncoeff*2*nexp+s[i]->jchip*np+2, ncoeff*2*nexp+s[i]->jchip*np+2) += Dx * Dx * isx2 + Dy * Dy * isy2;
	    }

	    // chip x star
	    a_data(ncoeff*2*nexp+s[i]->jchip*np  , size0+s[i]->jstar*2  ) -= Bx * s[i]->xi_a * isx2 + By * s[i]->eta_a * isy2;
	    a_data(ncoeff*2*nexp+s[i]->jchip*np  , size0+s[i]->jstar*2+1) -= Bx * s[i]->xi_d * isx2 + By * s[i]->eta_d * isy2;
	    a_data(ncoeff*2*nexp+s[i]->jchip*np+1, size0+s[i]->jstar*2  ) -= Cx * s[i]->xi_a * isx2 + Cy * s[i]->eta_a * isy2;
	    a_data(ncoeff*2*nexp+s[i]->jchip*np+1, size0+s[i]->jstar*2+1) -= Cx * s[i]->xi_d * isx2 + Cy * s[i]->eta_d * isy2;
	    a_data(size0+s[i]->jstar*2  , ncoeff*2*nexp+s[i]->jchip*np  ) -= Bx * s[i]->xi_a * isx2 + By * s[i]->eta_a * isy2;
	    a_data(size0+s[i]->jstar*2+1, ncoeff*2*nexp+s[i]->jchip*np  ) -= Bx * s[i]->xi_d * isx2 + By * s[i]->eta_d * isy2;
	    a_data(size0+s[i]->jstar*2  , ncoeff*2*nexp+s[i]->jchip*np+1) -= Cx * s[i]->xi_a * isx2 + Cy * s[i]->eta_a * isy2;
	    a_data(size0+s[i]->jstar*2+1, ncoeff*2*nexp+s[i]->jchip*np+1) -= Cx * s[i]->xi_d * isx2 + Cy * s[i]->eta_d * isy2;
	    if (allowRotation) {
		a_data(ncoeff*2*nexp+s[i]->jchip*np+2, size0+s[i]->jstar*2  ) -= Dx * s[i]->xi_a * isx2 + Dy * s[i]->eta_a * isy2;
		a_data(ncoeff*2*nexp+s[i]->jchip*np+2, size0+s[i]->jstar*2+1) -= Dx * s[i]->xi_d * isx2 + Dy * s[i]->eta_d * isy2;
		a_data(size0+s[i]->jstar*2  , ncoeff*2*nexp+s[i]->jchip*np+2) -= Dx * s[i]->xi_a * isx2 + Dy * s[i]->eta_a * isy2;
		a_data(size0+s[i]->jstar*2+1, ncoeff*2*nexp+s[i]->jchip*np+2) -= Dx * s[i]->xi_d * isx2 + Dy * s[i]->eta_d * isy2;
	    }

	    // star x star
	    a_data(size0+s[i]->jstar*2  , size0+s[i]->jstar*2  ) += s[i]->xi_a * s[i]->xi_a * isx2 + s[i]->eta_a * s[i]->eta_a * isy2;
	    a_data(size0+s[i]->jstar*2  , size0+s[i]->jstar*2+1) += s[i]->xi_a * s[i]->xi_d * isx2 + s[i]->eta_a * s[i]->eta_d * isy2;
	    a_data(size0+s[i]->jstar*2+1, size0+s[i]->jstar*2  ) += s[i]->xi_d * s[i]->xi_a * isx2 + s[i]->eta_d * s[i]->eta_a * isy2;
	    a_data(size0+s[i]->jstar*2+1, size0+s[i]->jstar*2+1) += s[i]->xi_d * s[i]->xi_d * isx2 + s[i]->eta_d * s[i]->eta_d * isy2;

	    b_data(ncoeff*2*nexp+s[i]->jchip*np  ) += Ax * Bx * isx2 + Ay * By * isy2;
	    b_data(ncoeff*2*nexp+s[i]->jchip*np+1) += Ax * Cx * isx2 + Ay * Cy * isy2;
	    if (allowRotation) {
		b_data(ncoeff*2*nexp+s[i]->jchip*np+2) += Ax * Dx * isx2 + Ay * Dy * isy2;
	    }

	    b_data(size0+2*s[i]->jstar  ) -= Ax * s[i]->xi_a * isx2 + Ay * s[i]->eta_a * isy2;
	    b_data(size0+2*s[i]->jstar+1) -= Ax * s[i]->xi_d * isx2 + Ay * s[i]->eta_d * isy2;
	}

	if (allowRotation) {
	    // \Sum d_theta = 0.0
	    for (int i = 0; i < nchip; i++) {
		a_data(ncoeff*2*nexp+i*np+2, ncoeff*2*nexp+nchip*np) = 1;
		a_data(ncoeff*2*nexp+nchip*np, ncoeff*2*nexp+i*np+2) = 1;
	    }
	}
    } else {
	for (int i = 0; i < nobs; i++) {
	    if (!o[i]->good) continue;
            ++numObsGood;
	    double Ax = o[i]->xi;
	    double Ay = o[i]->eta;
	    double Bx = 0.0;
	    double By = 0.0;
	    double Cx = 0.0;
	    double Cy = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		pu(k) = pow(o[i]->u, xorder[k]);
		pv(k) = pow(o[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu(k)   * pv(k);
		Ay -= b[o[i]->iexp][k] * pu(k)   * pv(k);
		Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		Cx += a[o[i]->iexp][k] * pu(k)   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[o[i]->iexp][k] * pu(k)   * pow(o[i]->v, yorder[k]-1) * yorder[k];
	    }
	    double dxi  = Bx * o[i]->xerr + Cx * o[i]->yerr;
	    double deta = By * o[i]->xerr + Cy * o[i]->yerr;
	    isx2 = 1.0 / (pow(dxi,  2) + pow(catRMS, 2));
	    isy2 = 1.0 / (pow(deta, 2) + pow(catRMS, 2));

	    for (int k = 0; k < ncoeff; k++) {
		b_data(k+       ncoeff*2*o[i]->jexp) += Ax * pu(k) * pv(k) * isx2;
		b_data(k+ncoeff+ncoeff*2*o[i]->jexp) += Ay * pu(k) * pv(k) * isy2;
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data(j+       ncoeff*2*o[i]->jexp, k+       ncoeff*2*o[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isx2;
		    a_data(j+ncoeff+ncoeff*2*o[i]->jexp, k+ncoeff+ncoeff*2*o[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isy2;
		}
	    }
	}

	for (int i = 0; i < nSobs; i++) {
	    if (!s[i]->good || s[i]->jstar == -1) continue;
            ++numStarGood;
	    double Ax = s[i]->xi;
	    double Ay = s[i]->eta;
	    double Bx = 0.0;
	    double By = 0.0;
	    double Cx = 0.0;
	    double Cy = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		pu(k) = pow(s[i]->u, xorder[k]);
		pv(k) = pow(s[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[s[i]->iexp][k] * pu(k)   * pv(k);
		Ay -= b[s[i]->iexp][k] * pu(k)   * pv(k);
		Bx += a[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		By += b[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pv(k)   * xorder[k];
		Cx += a[s[i]->iexp][k] * pu(k)   * pow(s[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[s[i]->iexp][k] * pu(k)   * pow(s[i]->v, yorder[k]-1) * yorder[k];
	    }
	    double dxi  = Bx * s[i]->xerr + Cx * s[i]->yerr;
	    double deta = By * s[i]->xerr + Cy * s[i]->yerr;
	    isx2 = 1.0 / pow(dxi,  2);
	    isy2 = 1.0 / pow(deta, 2);

	    for (int k = 0; k < ncoeff; k++) {
		b_data(k+       ncoeff*2*s[i]->jexp) += Ax * pu(k) * pv(k) * isx2;
		b_data(k+ncoeff+ncoeff*2*s[i]->jexp) += Ay * pu(k) * pv(k) * isy2;
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data(j+       ncoeff*2*s[i]->jexp, k+       ncoeff*2*s[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isx2;
		    a_data(j+ncoeff+ncoeff*2*s[i]->jexp, k+ncoeff+ncoeff*2*s[i]->jexp) += pu(j) * pv(j) * pu(k) * pv(k) * isy2;
		}

		// coeff x star
		a_data(k+       ncoeff*2*s[i]->jexp, size0+s[i]->jstar*2  ) -= s[i]->xi_a  * pu(k) * pv(k) * isx2;
		a_data(k+       ncoeff*2*s[i]->jexp, size0+s[i]->jstar*2+1) -= s[i]->xi_d  * pu(k) * pv(k) * isx2;
		a_data(k+ncoeff+ncoeff*2*s[i]->jexp, size0+s[i]->jstar*2  ) -= s[i]->eta_a * pu(k) * pv(k) * isy2;
		a_data(k+ncoeff+ncoeff*2*s[i]->jexp, size0+s[i]->jstar*2+1) -= s[i]->eta_d * pu(k) * pv(k) * isy2;
		a_data(size0+s[i]->jstar*2  , k+       ncoeff*2*s[i]->jexp) -= s[i]->xi_a  * pu(k) * pv(k) * isx2;
		a_data(size0+s[i]->jstar*2+1, k+       ncoeff*2*s[i]->jexp) -= s[i]->xi_d  * pu(k) * pv(k) * isx2;
		a_data(size0+s[i]->jstar*2  , k+ncoeff+ncoeff*2*s[i]->jexp) -= s[i]->eta_a * pu(k) * pv(k) * isy2;
		a_data(size0+s[i]->jstar*2+1, k+ncoeff+ncoeff*2*s[i]->jexp) -= s[i]->eta_d * pu(k) * pv(k) * isy2;
	    }

	    // star x star
	    a_data(size0+s[i]->jstar*2  , size0+s[i]->jstar*2  ) += s[i]->xi_a * s[i]->xi_a * isx2 + s[i]->eta_a * s[i]->eta_a * isy2;
	    a_data(size0+s[i]->jstar*2  , size0+s[i]->jstar*2+1) += s[i]->xi_a * s[i]->xi_d * isx2 + s[i]->eta_a * s[i]->eta_d * isy2;
	    a_data(size0+s[i]->jstar*2+1, size0+s[i]->jstar*2  ) += s[i]->xi_d * s[i]->xi_a * isx2 + s[i]->eta_d * s[i]->eta_a * isy2;
	    a_data(size0+s[i]->jstar*2+1, size0+s[i]->jstar*2+1) += s[i]->xi_d * s[i]->xi_d * isx2 + s[i]->eta_d * s[i]->eta_d * isy2;

	    b_data(size0+2*s[i]->jstar  ) -= Ax * s[i]->xi_a * isx2 + Ay * s[i]->eta_a * isy2;
	    b_data(size0+2*s[i]->jstar+1) -= Ax * s[i]->xi_d * isx2 + Ay * s[i]->eta_d * isy2;
	}
    }

    std::cout << "Number good: " << numObsGood << ", " << numStarGood << std::endl;

    Eigen::VectorXd coeff = solveMatrix(size, a_data, b_data);
//    Eigen::VectorXd coeff = a_data.partialPivLu().solve(b_data);

    return coeff;
}

double calcChi2(std::vector<Obs::Ptr>& o, Coeff::Ptr c, Poly::Ptr p)
{
    int nobs  = o.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double *a = c->a;
    double *b = c->b;

    double chi2 = 0.0;
    for (int i = 0; i < nobs; i++) {
	if (!o[i]->good) continue;
	double Ax = o[i]->xi;
	double Ay = o[i]->eta;
	for (int k = 0; k < ncoeff; k++) {
	    Ax -= a[k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Ay -= b[k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	}
	chi2 += Ax * Ax + Ay * Ay;
    }

    return chi2;
}

double calcChi2(std::vector<Obs::Ptr>& o, CoeffSet& coeffVec, Poly::Ptr p, bool norm=false)
{
    int nobs  = o.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    std::map<int, double*> a;
    std::map<int, double*> b;
    for (CoeffSet::iterator it = coeffVec.begin(); it != coeffVec.end(); it++) {
	a.insert(std::map<int, double*>::value_type(it->first, it->second->a));
	b.insert(std::map<int, double*>::value_type(it->first, it->second->b));
    }

    double chi2 = 0.0;
    int num = 0;
    for (int i = 0; i < nobs; i++) {
	if (!o[i]->good) continue;
	double Ax = o[i]->xi;
	double Ay = o[i]->eta;
	for (int k = 0; k < ncoeff; k++) {
	    Ax -= a[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Ay -= b[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	}
	chi2 += Ax * Ax + Ay * Ay;
	num++;
    }

    if (norm)
	return chi2/num;
    else
	return chi2;
}

void flagObj2(std::vector<Obs::Ptr>& o, CoeffSet& coeffVec, Poly::Ptr p, double e2, double catRMS=0.0)
{
    int nobs  = o.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    std::map<int, double*> a;
    std::map<int, double*> b;
    for (CoeffSet::iterator it = coeffVec.begin(); it != coeffVec.end(); it++) {
	a.insert(std::map<int, double*>::value_type(it->first, it->second->a));
	b.insert(std::map<int, double*>::value_type(it->first, it->second->b));
    }

    int nreject = 0;
    for (int i = 0; i < nobs; i++) {
	if (!o[i]->good) continue;
	double Ax = o[i]->xi;
	double Ay = o[i]->eta;
	double Bx = 0.0;
	double By = 0.0;
	double Cx = 0.0;
	double Cy = 0.0;
	for (int k = 0; k < ncoeff; k++) {
	    Ax -= a[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Ay -= b[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]) * xorder[k];
	    By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]) * xorder[k];
	    Cx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]-1) * yorder[k];
	    Cy += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]-1) * yorder[k];
	}
	//double dxi  = Bx * o[i]->xerr + Cx * o[i]->yerr;
	//double deta = By * o[i]->xerr + Cy * o[i]->yerr;
	//double chi2 = Ax * Ax / (dxi * dxi + catRMS * catRMS) + Ay * Ay / (deta * deta + catRMS * catRMS);
	//if (chi2 > e2) {
	double chi2 = Ax * Ax + Ay * Ay;
	if (chi2 > e2+9.0*catRMS*catRMS) {
	    o[i]->good = false;
	    nreject++;
	} else {
	    //o[i]->good = true;
	}
    }
    printf("nreject = %d\n", nreject);

}

double calcChi2_Star(std::vector<Obs::Ptr>& o, std::vector<Obs::Ptr>& s, CoeffSet& coeffVec, Poly::Ptr p)
{
    double chi2 = 0.0;
    chi2 += calcChi2(o, coeffVec, p);
    chi2 += calcChi2(s, coeffVec, p);

    return chi2;
}

int lsst::meas::mosaic::flagSuspect(SourceGroup &allMat,
				    SourceGroup &allSource,
				    WcsDic &wcsDic)
{
    std::vector<int> visits;
    for (WcsDic::iterator it = wcsDic.begin(); it != wcsDic.end(); it++) {
	visits.push_back(it->first);
    }

    for (size_t j = 0; j < visits.size()-1; j++) {
	int visit_ref = visits[j];
	for (size_t i = j+1; i < visits.size(); i++) {
	    int visit_targ = visits[i];
	    std::vector<PTR(Source)> refs;
	    std::vector<PTR(Source)> targs;
	    std::vector<double> mref;
	    std::vector<double> mtarg;
	    for (size_t l = 0; l < allMat.size(); l++) {
		std::vector<PTR(Source)> mm = allMat[l];
		int j_ref = -1;
		int j_targ = -1;
		for (size_t k = 1; k < mm.size(); k++) {
		    if (mm[k]->getExp() == visit_ref) {
			j_ref = k;
		    } else if (mm[k]->getExp() == visit_targ) {
			j_targ = k;
		    }
		}
		if (j_ref != -1 && j_targ != -1 &&
		    mm[j_ref]->getFlux() > 0.0 &&
		    mm[j_targ]->getFlux() > 0.0) {
		    refs.push_back(mm[j_ref]);
		    targs.push_back(mm[j_targ]);
		    mref.push_back(-2.5*log10(mm[j_ref]->getFlux()));
		    mtarg.push_back(-2.5*log10(mm[j_targ]->getFlux()));
		}
	    }
	    for (size_t l = 0; l < allSource.size(); l++) {
		std::vector<PTR(Source)> ss = allSource[l];
		int j_ref = -1;
		int j_targ = -1;
		for (size_t k = 1; k < ss.size(); k++) {
		    if (ss[k]->getExp() == visit_ref) {
			j_ref = k;
		    } else if (ss[k]->getExp() == visit_targ) {
			j_targ = k;
		    }
		}
		if (j_ref != -1 && j_targ != -1 &&
		    ss[j_ref]->getFlux() > 0.0 &&
		    ss[j_targ]->getFlux() > 0.0) {
		    refs.push_back(ss[j_ref]);
		    targs.push_back(ss[j_targ]);
		    mref.push_back(-2.5*log10(ss[j_ref]->getFlux()));
		    mtarg.push_back(-2.5*log10(ss[j_targ]->getFlux()));
		}
	    }

	    if (mref.size() < 10) {
		printf("%d %d\n", visit_ref, visit_targ);
		continue;
	    }

	    std::vector<double> dm;
	    for (size_t k = 0; k < mref.size(); k++) {
		dm.push_back(mtarg[k] - mref[k]);
	    }
	    std::sort(dm.begin(), dm.end());
	    double med = dm[dm.size()/2];
	    double Q1 = dm[dm.size()/10];
	    double Q3 = dm[dm.size()*9/10];
	    double SIQR = 0.5 * (Q3 - Q1);

	    int ngood = 0;
	    int nbad = 0;
	    for (size_t k = 0; k < mref.size(); k++) {
		if (fabs(mtarg[k]-mref[k]-med) > 3.0 * SIQR) {
		    refs[k]->setFlux(-9999);
		    targs[k]->setFlux(-9999);
		    nbad += 1;
		} else {
		    ngood += 1;
		}
	    }
	    printf("%d %d %6.3f %5.3f %5d %5d\n", visit_ref, visit_targ, med, SIQR, ngood, nbad);
	}
    }

    return 0;
}

ObsVec
lsst::meas::mosaic::obsVecFromSourceGroup(SourceGroup const &all,
					  WcsDic &wcsDic,
					  CcdSet &ccdSet)
{
    std::vector<Obs::Ptr> obsVec;
    for (size_t i = 0; i < all.size(); i++) {
        SourceSet ss = all[i];
	double ra  = ss[0]->getRa().asRadians();
	double dec = ss[0]->getDec().asRadians();
	double mag_cat;
	double err_cat;
	if (ss[0]->getFlux() > 0.0 && ss[0]->getFluxErr() > 0.0) {
	    mag_cat = -2.5*log10(ss[0]->getFlux());
	    err_cat = 2.5/M_LN10*ss[0]->getFluxErr()/ss[0]->getFlux();
	} else {
	    mag_cat = -9999;
	    err_cat = -9999;
	}
	//std::cout << ra << " " << dec << std::endl;
	//std::cout << mag0 << std::endl;
	for (size_t j = 1; j < ss.size(); j++) {
	    long id   = ss[j]->getId();
	    int iexp  = ss[j]->getExp();
	    int ichip = ss[j]->getChip();
	    double x = ss[j]->getX();
	    double y = ss[j]->getY();
	    Obs::Ptr o = Obs::Ptr(new Obs(id, ra, dec, x, y, ichip, iexp));

	    int jexp = 0;
	    for (WcsDic::iterator it = wcsDic.begin(); it != wcsDic.end(); it++) {
		if (it->first == iexp) break;
		jexp++;
	    }
	    o->jexp = jexp;

	    int jchip = 0;
	    for (CcdSet::iterator it = ccdSet.begin(); it != ccdSet.end(); it++) {
		if (it->first == ichip) break;
		jchip++;
	    }
	    o->jchip = jchip;

	    o->mag_cat = mag_cat;
	    o->err_cat = err_cat;
	    o->mag0 = mag_cat;
	    lsst::afw::geom::PointD crval
		= wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::geom::radians);
	    o->setXiEta(crval[0], crval[1]);
	    o->setUV(ccdSet[ichip]);
	    o->xerr = ss[j]->getXErr();
	    o->yerr = ss[j]->getYErr();
	    if (lsst::utils::isnan(o->xerr) || lsst::utils::isnan(o->yerr))
		o->good = false;
	    //printf("%d %e %e\n", o->good, o->xerr, o->yerr);
	    o->istar = i;
	    if (ss[0]->getAstromBad() || ss[j]->getAstromBad()) {
		o->good = false;
	    }
	    if (ss[j]->getFlux() > 0.0 && ss[j]->getFluxErr() > 0.0 && ss[j]->getFlux() < 1.0E+10) {
		o->mag = -2.5*log10(ss[j]->getFlux());
		o->err = 2.5/M_LN10*ss[j]->getFluxErr()/ss[j]->getFlux();
	    } else {
		o->mag = -9999;
		o->err = -9999;
		o->good = false;
	    }
	    obsVec.push_back(o);
	}
    }

    return obsVec;
}

Eigen::VectorXd solveSIP_P(Poly::Ptr p,
			   std::vector<Obs::Ptr> &obsVec) {
    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    Eigen::MatrixXd a_data = Eigen::MatrixXd::Zero(ncoeff, ncoeff);
    Eigen::MatrixXd d_data = Eigen::MatrixXd::Zero(ncoeff, ncoeff);
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(ncoeff);
    Eigen::VectorXd c_data = Eigen::VectorXd::Zero(ncoeff);

    Eigen::VectorXd pu(ncoeff);
    Eigen::VectorXd pv(ncoeff);

    for (size_t k = 0; k < obsVec.size(); k++) {
	Obs::Ptr o = obsVec[k];
	if (o->good) {
	    for (int j = 0; j < ncoeff; j++) {
		pu(j) = pow(o->U, xorder[j]);
		pv(j) = pow(o->V, yorder[j]);
	    }
	    for (int j = 0; j < ncoeff; j++) {
		b_data(j) += (o->u - o->U) * pu(j) * pv(j);
		c_data(j) += (o->v - o->V) * pu(j) * pv(j);
		for (int i = 0; i < ncoeff; i++) {
		    a_data(i, j) += pu(j) * pv(j) * pu(i) * pv(i);
		    d_data(i, j) += pu(j) * pv(j) * pu(i) * pv(i);
		}
	    }
	}
    }

    Eigen::VectorXd coeffA = solveMatrix(ncoeff, a_data, b_data);
    Eigen::VectorXd coeffB = solveMatrix(ncoeff, d_data, c_data);

    Eigen::VectorXd coeff(2*ncoeff);
    for (int i = 0; i < ncoeff; i++) {
	coeff(i)        = coeffA(i);
	coeff(i+ncoeff) = coeffB(i);
    }

    return coeff;
}

void setCRVALtoDetJPeak(Coeff::Ptr c) {
    double w = (3.-sqrt(5.))/2.;
    double ua, ub, uc, ux;
    double va, vb, vc, vx;
    double fb, fx;
    double u, v, upre, vpre;

    u = upre = 0.0;
    v = vpre = 0.0;

    for (int i = 0; i < 10; i++) {
	ua = u - 3000 / pow(2, i);
	uc = u + 3000 / pow(2, i);
	ub = ua * (1-w) + uc * w;

	//fa = c->detJ(ua, v);
	fb = c->detJ(ub, v);
	//fc = c->detJ(uc, v);

	while (1) {
	    if (uc - ub > ub - ua) {
		ux = ub * (1-w) + uc * w;
	    } else {
		ux = ua * (1-w) + ub * w;
	    }
	    fx = c->detJ(ux, v);
	    if (uc - ub > ub - ua) {
		if (fx > fb) {
		    ua = ub;
		    ub = ux;
		    //fa = c->detJ(ua, v);
		    fb = c->detJ(ub, v);
		} else {
		    uc = ux;
		    //fc = c->detJ(uc, v);
		}
	    } else {
		if (fx > fb) {
		    uc = ub;
		    ub = ux;
		    //fc = c->detJ(uc, v);
		    fb = c->detJ(ub, v);
		} else {
		    ua = ux;
		    //fa = c->detJ(ua, v);
		}
	    }
	    if (uc - ua < 0.01) break;
	}

	u = ub;

	va = v - 3000 / pow(2, i);
	vc = v + 3000 / pow(2, i);
	vb = va * (1-w) + vc * w;

	//fa = c->detJ(u, va);
	fb = c->detJ(u, vb);
	//fc = c->detJ(u, vc);

	while (1) {
	    if (vc - vb > vb - va) {
		vx = vb * (1-w) + vc * w;
	    } else {
		vx = va * (1-w) + vb * w;
	    }
	    fx = c->detJ(u, vx);
	    if (vc - vb > vb - va) {
		if (fx > fb) {
		    va = vb;
		    vb = vx;
		    //fa = c->detJ(u, va);
		    fb = c->detJ(u, vb);
		} else {
		    vc = vx;
		    //fc = c->detJ(u, vc);
		}
	    } else {
		if (fx > fb) {
		    vc = vb;
		    vb = vx;
		    //fc = c->detJ(u, vc);
		    fb = c->detJ(u, vb);
		} else {
		    va = vx;
		    //fa = c->detJ(u, va);
		}
	    }
	    if (vc - va < 0.01) break;
	}

	v = vb;

	if (fabs(u-upre) < 0.01 && fabs(v-vpre) < 0.01) break;
    }

    double xi, eta;
    c->uvToXiEta(u, v, &xi, &eta);
    xi  = xi  * D2R;
    eta = eta * D2R;

    double phi, theta;
    phi = atan2(xi, eta);
    theta = atan2(1.0, sqrt(xi*xi+eta*eta));

    double x = sin(theta);
    double y = cos(theta)*sin(phi);
    double z = cos(theta)*cos(phi);

    double alpha = atan2(y, z*sin(c->D)-x*cos(c->D));
    if (z*sin(c->D)-x*cos(c->D) < 0.0) alpha += M_PI;
    if (alpha > M_PI) alpha -= 2*M_PI;
    double sinalpha = sin(alpha);
    double delta = atan2(x*sin(c->D)+z*cos(c->D), -y/sinalpha);
    alpha = -alpha + c->A;
    if (alpha > 2*M_PI) alpha -= 2*M_PI;

    //printf("%f %f\n", c->A*R2D, c->D*R2D);
    //printf("%f %f\n", alpha*R2D, delta*R2D);

    c->A = alpha;
    c->D = delta;
}

CoeffSet
initialFit(int nexp,
	   ObsVec &matchVec,
	   WcsDic &wcsDic,
	   CcdSet &ccdSet,
	   Poly::Ptr &p) {
    int nMobs = matchVec.size();

    // Solve for polynomial coefficients and crvals
    // for each exposure separately
    // These values will be used as initial guess for
    // the subsequent fitting

    CoeffSet coeffVec;

    for (WcsDic::iterator it =  wcsDic.begin(); it != wcsDic.end(); it++) {
	int iexp = it->first;

	// Select objects for a specific exposure id
	std::vector<Obs::Ptr> obsVec_sub;
	for (int j = 0; j < nMobs; j++) {
	    if (matchVec[j]->iexp == iexp) {
		obsVec_sub.push_back(matchVec[j]);
	    }
	}

	// Solve for polinomial and crval
	Eigen::VectorXd a = solveForCoeff(obsVec_sub, p);

	double chi2 = calcChi(obsVec_sub, a, p);
	printf("visit: %d\n", iexp);
	printf("calcChi: %e\n", chi2);
	double e2 = chi2 / obsVec_sub.size();
	flagObj(obsVec_sub, a, p, 9.0*e2);

	a = solveForCoeff(obsVec_sub, p);
	chi2 = calcChi(obsVec_sub, a, p);
	printf("calcChi: %e\n", chi2);

	// Store solution into Coeff class
	Coeff::Ptr c = Coeff::Ptr(new Coeff(p));
	c->iexp = iexp;
	for (int k = 0; k < p->ncoeff; k++) {
	    c->a[k] = a(k);
	    c->b[k] = a(k+p->ncoeff);
	}
	lsst::afw::geom::PointD crval
	    = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::geom::radians);
	c->A = crval[0] + a(p->ncoeff*2);
	c->D = crval[1] + a(p->ncoeff*2+1);
	c->x0 = c->y0 = 0.0;

	for (size_t j = 0; j < obsVec_sub.size(); j++) {
	    obsVec_sub[j]->setXiEta(c->A, c->D);
	}

	a = solveForCoeffWithOffset(obsVec_sub, c, p);

	// Store solution into Coeff class
	for (int k = 0; k < p->ncoeff; k++) {
	    c->a[k] += a(k);
	    c->b[k] += a(k+p->ncoeff);
	}
	c->x0 += a(2*p->ncoeff);
	c->y0 += a(2*p->ncoeff+1);

	for (size_t j = 0; j < obsVec_sub.size(); j++) {
	    obsVec_sub[j]->setUV(ccdSet[obsVec_sub[j]->ichip], c->x0, c->y0);
	}
	chi2 = calcChi2(obsVec_sub, c, p);
	printf("calcChi2: %e\n", chi2);

//	setCRVALtoDetJPeak(c);

	for (size_t j = 0; j < obsVec_sub.size(); j++) {
	    obsVec_sub[j]->setXiEta(c->A, c->D);
	}

	a = solveForCoeffWithOffset(obsVec_sub, c, p);

	// Store solution into Coeff class
	for (int k = 0; k < p->ncoeff; k++) {
	    c->a[k] += a(k);
	    c->b[k] += a(k+p->ncoeff);
	}
	c->x0 += a(2*p->ncoeff);
	c->y0 += a(2*p->ncoeff+1);

	for (size_t j = 0; j < obsVec_sub.size(); j++) {
	    obsVec_sub[j]->setUV(ccdSet[obsVec_sub[j]->ichip], c->x0, c->y0);
	}
	chi2 = calcChi2(obsVec_sub, c, p);
	printf("calcChi2: %e\n", chi2);

	/////////////////////////////////////////////////////////////////////////////////
	a = solveForCoeffWithOffset(obsVec_sub, c, p);

	// Store solution into Coeff class
	for (int k = 0; k < p->ncoeff; k++) {
	    c->a[k] += a(k);
	    c->b[k] += a(k+p->ncoeff);
	}
	c->x0 += a(2*p->ncoeff);
	c->y0 += a(2*p->ncoeff+1);

	for (size_t j = 0; j < obsVec_sub.size(); j++) {
	    obsVec_sub[j]->setUV(ccdSet[obsVec_sub[j]->ichip], c->x0, c->y0);
	}
	chi2 = calcChi2(obsVec_sub, c, p);
	printf("calcChi2: %e\n", chi2);
	/////////////////////////////////////////////////////////////////////////////////

	coeffVec.insert(std::map<int, Coeff::Ptr>::value_type(iexp, c));

    }

    return coeffVec;
}

CoeffSet
lsst::meas::mosaic::solveMosaic_CCD_shot(int order,
					 int nmatch,
					 ObsVec &matchVec,
					 WcsDic &wcsDic,
					 CcdSet &ccdSet,
					 bool solveCcd,
					 bool allowRotation,
					 bool verbose,
					 double catRMS,
					 bool writeSnapshots,
					 std::string const & snapshotDir
)
{
    boost::filesystem::path snapshotPath(snapshotDir);

    Poly::Ptr p = Poly::Ptr(new Poly(order));

    int nMobs = matchVec.size();

    int nexp = wcsDic.size();
    int nchip = ccdSet.size();
    int ncoeff = p->ncoeff;

    if (writeSnapshots) {
        writeObsVec((snapshotPath / "match-initial-0.fits").native(), matchVec);
    }

    // Solve for polynomial coefficients and crvals
    // for each exposure separately
    // These values will be used as initial guess for
    // the subsequent fitting

    CoeffSet coeffVec = initialFit(nexp, matchVec, wcsDic, ccdSet, p);

    // Update Xi and Eta using new crval (rac and decc)
    for (int i = 0; i < nMobs; i++) {
	double rac  = coeffVec[matchVec[i]->iexp]->A;
	double decc = coeffVec[matchVec[i]->iexp]->D;
	matchVec[i]->setXiEta(rac, decc);
	matchVec[i]->setFitVal(coeffVec[matchVec[i]->iexp], p);
    }

    if (writeSnapshots) {
        writeObsVec((snapshotPath / "match-initial-1.fits").native(), matchVec);
    }

    printf("Before fitting calcChi2: %e\n",
	   calcChi2(matchVec, coeffVec, p));
    printf("Before fitting matched: %5.3f (arcsec)\n",
	   sqrt(calcChi2(matchVec, coeffVec, p, true))*3600.0);

    for (int k = 0; k < 3; k++) {
	Eigen::VectorXd coeff = solveLinApprox(matchVec, coeffVec, nchip, p, solveCcd, allowRotation, catRMS);

	int j = 0;
	for (CoeffSet::iterator it = coeffVec.begin(); it != coeffVec.end(); it++, j++) {
	    for (int i = 0; i < ncoeff; i++) {
		it->second->a[i] += coeff(2*ncoeff*j+i);
		it->second->b[i] += coeff(2*ncoeff*j+i+ncoeff);
	    }
	}

	if (solveCcd) {
	    if (allowRotation) {
		int i = 0;
		for (CcdSet::iterator it = ccdSet.begin(); it != ccdSet.end(); it++, i++) {
            // offset is calculated in pixel coordinates.
            afw::geom::Extent2D offset(coeff(2*ncoeff*nexp+3*i), coeff(2*ncoeff*nexp+3*i+1));

            // but must be transformed to mm to be applied.
            auto pixelSize = it->second->getPixelSize();
            auto scaling = afw::geom::LinearTransform::makeScaling(pixelSize.getX(), pixelSize.getY());
            offset = scaling(offset);

            afw::cameraGeom::Orientation newOrientation(
                it->second->getOrientation().getFpPosition() + offset,
                it->second->getOrientation().getReferencePoint(),
                it->second->getOrientation().getYaw() + coeff(2*ncoeff*nexp + 3*i + 2) * afw::geom::radians,
                it->second->getOrientation().getPitch(),
                it->second->getOrientation().getRoll());

            afw::cameraGeom::CameraTransformMap::Transforms newTr;

            // Transform from pixels to pixels is the same on the new Detector as the old.
            newTr[afw::cameraGeom::CameraSys(afw::cameraGeom::PIXELS, it->second->getName())] =
                 it->second->getTransformMap()[afw::cameraGeom::CameraSys(afw::cameraGeom::PIXELS,
                                                                          it->second->getName())];

            // Transform from pixels to focal plane has to be recalculated.
            newTr[afw::cameraGeom::FOCAL_PLANE] = boost::make_shared<afw::geom::AffineXYTransform const>
                                                  (newOrientation.makePixelFpTransform(pixelSize));

            // We should not require any other transformations within meas_mosaic.

            // Replace the old Detector with a shifted one.
            it->second = boost::make_shared<afw::cameraGeom::Detector>(
                it->second->getName(),
                it->second->getId(),
                it->second->getType(),
                it->second->getSerial(),
                it->second->getBBox(),
                it->second->getAmpInfoCatalog(),
                newOrientation,
                it->second->getPixelSize(),
                newTr);

		}
	    } else {
		int i = 0;
		for (CcdSet::iterator it = ccdSet.begin(); it != ccdSet.end(); it++, i++) {
            // Logic as above, but with no change to Detector yaw.
            afw::geom::Extent2D offset(coeff(2*ncoeff*nexp+2*i), coeff(2*ncoeff*nexp+2*i+1));
            auto pixelSize = it->second->getPixelSize();
            auto scaling = afw::geom::LinearTransform::makeScaling(pixelSize.getX(), pixelSize.getY());
            offset = scaling(offset);

            afw::cameraGeom::Orientation newOrientation(
                it->second->getOrientation().getFpPosition() + offset,
                it->second->getOrientation().getReferencePoint(),
                it->second->getOrientation().getYaw(),
                it->second->getOrientation().getPitch(),
                it->second->getOrientation().getRoll());

            afw::cameraGeom::CameraTransformMap::Transforms newTr;

            // Transform from pixels to pixels is the same on the new Detector as the old.
            newTr[afw::cameraGeom::CameraSys(afw::cameraGeom::PIXELS, it->second->getName())] =
                 it->second->getTransformMap()[afw::cameraGeom::CameraSys(afw::cameraGeom::PIXELS,
                                                                          it->second->getName())];

            // Transform from pixels to focal plane has to be recalculated.
            newTr[afw::cameraGeom::FOCAL_PLANE] = boost::make_shared<afw::geom::AffineXYTransform const>
                                                  (newOrientation.makePixelFpTransform(pixelSize));

            // We should not require any other transformations within meas_mosaic.

            // Replace the old Detector with a shifted one.
            it->second = boost::make_shared<afw::cameraGeom::Detector>(
                it->second->getName(),
                it->second->getId(),
                it->second->getType(),
                it->second->getSerial(),
                it->second->getBBox(),
                it->second->getAmpInfoCatalog(),
                newOrientation,
                it->second->getPixelSize(),
                newTr);
		}
	    }
	}

	for (int i = 0; i < nMobs; i++) {
	    matchVec[i]->setUV(ccdSet[matchVec[i]->ichip], coeffVec[matchVec[i]->iexp]->x0, coeffVec[matchVec[i]->iexp]->y0);
	    matchVec[i]->setFitVal(coeffVec[matchVec[i]->iexp], p);
	}

	if (writeSnapshots) {
	    writeObsVec((snapshotPath / (boost::format("match-iter-%d.fits") % k).str()).native(), matchVec);
	}

	printf("%dth iteration calcChi2: %e\n", (k+1), calcChi2(matchVec, coeffVec, p));
	printf("%dth iteration matched: %5.3f (arcsec)\n",
	       (k+1),
	       sqrt(calcChi2(matchVec, coeffVec, p, true))*3600.0);
	double drm = calcChi2(matchVec, coeffVec, p, true);
	flagObj2(matchVec, coeffVec, p, 9.0*drm, catRMS);
    }

    std::map<int, Eigen::Matrix2d> cd;
    for (CoeffSet::iterator it = coeffVec.begin(); it != coeffVec.end(); it++) {
	Eigen::Matrix2d c;
	c << it->second->a[0], it->second->a[1], it->second->b[0], it->second->b[1];
	cd.insert(std::map<int, Eigen::Matrix2d>::value_type(it->first, c));
    }
    for (int i = 0; i < nMobs; i++) {
	double CD1_1 = cd[matchVec[i]->iexp](0,0);
	double CD1_2 = cd[matchVec[i]->iexp](0,1);
	double CD2_1 = cd[matchVec[i]->iexp](1,0);
	double CD2_2 = cd[matchVec[i]->iexp](1,1);
	double det = CD1_1 * CD2_2 - CD1_2 * CD2_1;
	matchVec[i]->U = ( matchVec[i]->xi * CD2_2 - matchVec[i]->eta * CD1_2) / det;
	matchVec[i]->V = (-matchVec[i]->xi * CD2_1 + matchVec[i]->eta * CD1_1) / det;
    }

    for (CoeffSet::iterator it = coeffVec.begin(); it != coeffVec.end(); it++) {
	std::vector<Obs::Ptr> obsVec_sub;
	for (size_t j = 0; j < matchVec.size(); j++) {
	    Obs::Ptr iobs = matchVec[j];
	    if (iobs->iexp == it->first) {
		obsVec_sub.push_back(iobs);
	    }
	}
	Eigen::VectorXd a = solveSIP_P(p, obsVec_sub);
	for (int k = 0; k < p->ncoeff; k++) {
	    it->second->ap[k] = a(k);
	    it->second->bp[k] = a(k+p->ncoeff);
	}
    }

    for (int i = 0; i < nMobs; i++) {
	matchVec[i]->setFitVal2(coeffVec[matchVec[i]->iexp], p);
    }

    return coeffVec;
}

CoeffSet
lsst::meas::mosaic::solveMosaic_CCD(int order,
				    int nmatch,
				    int nsource,
				    ObsVec &matchVec,
				    ObsVec &sourceVec,
				    WcsDic &wcsDic,
				    CcdSet &ccdSet,
				    bool solveCcd,
				    bool allowRotation,
				    bool verbose,
				    double catRMS,
				    bool writeSnapshots,
				    std::string const & snapshotDir
)
{
    boost::filesystem::path snapshotPath(snapshotDir);

    Poly::Ptr p = Poly::Ptr(new Poly(order));

    int nMobs = matchVec.size();
    int nSobs = sourceVec.size();

    int nexp = wcsDic.size();
    int nchip = ccdSet.size();
    int ncoeff = p->ncoeff;
    int nstar = nsource;

    if (writeSnapshots) {
        writeObsVec((snapshotPath / "match-initial-0.fits").native(), matchVec);
        writeObsVec((snapshotPath / "source-initial-0.fits").native(), sourceVec);
    }

    // Solve for polynomial coefficients and crvals
    // for each exposure separately
    // These values will be used as initial guess for
    // the subsequent fitting

    CoeffSet coeffVec = initialFit(nexp, matchVec, wcsDic, ccdSet, p);

    // Update (xi, eta) and (u, v) using initial fitting resutls
    for (int i = 0; i < nMobs; i++) {
	double rac  = coeffVec[matchVec[i]->iexp]->A;
	double decc = coeffVec[matchVec[i]->iexp]->D;
	matchVec[i]->setXiEta(rac, decc);
	matchVec[i]->setUV(ccdSet[matchVec[i]->ichip],
			   coeffVec[matchVec[i]->iexp]->x0, coeffVec[matchVec[i]->iexp]->y0);
	matchVec[i]->setFitVal(coeffVec[matchVec[i]->iexp], p);
    }
    for (int i = 0; i < nSobs; i++) {
	double rac  = coeffVec[sourceVec[i]->iexp]->A;
	double decc = coeffVec[sourceVec[i]->iexp]->D;
	sourceVec[i]->setXiEta(rac, decc);
	sourceVec[i]->setUV(ccdSet[sourceVec[i]->ichip],
			    coeffVec[sourceVec[i]->iexp]->x0, coeffVec[sourceVec[i]->iexp]->y0);
	sourceVec[i]->setFitVal(coeffVec[sourceVec[i]->iexp], p);
    }

    if (writeSnapshots) {
        writeObsVec((snapshotPath / "match-initial-1.fits").native(), matchVec);
        writeObsVec((snapshotPath / "source-initial-1.fits").native(), sourceVec);
    }

    printf("Before fitting calcChi2: %e %e\n",
	   calcChi2(matchVec, coeffVec, p),
	   calcChi2_Star(matchVec, sourceVec, coeffVec, p));
    printf("Before fitting matched: %5.3f (arcsec) sources: %5.3f (arcsec)\n",
	   sqrt(calcChi2(matchVec, coeffVec, p, true))*3600.0,
	   sqrt(calcChi2(sourceVec, coeffVec, p, true))*3600.0);

    for (int k = 0; k < 3; k++) {
	Eigen::VectorXd coeff = solveLinApprox_Star(matchVec, sourceVec, nstar, coeffVec, nchip, p, solveCcd, allowRotation, catRMS);

	int j = 0;
	for (CoeffSet::iterator it = coeffVec.begin(); it != coeffVec.end(); it++, j++) {
	    for (int i = 0; i < ncoeff; i++) {
		it->second->a[i] += coeff(2*ncoeff*j+i);
		it->second->b[i] += coeff(2*ncoeff*j+i+ncoeff);
	    }
	}

	if (solveCcd) {
	    if (allowRotation) {
		int i = 0;
		for (CcdSet::iterator it = ccdSet.begin(); it != ccdSet.end(); it++, i++) {
            // offset is calculated in pixel coordinates.
            afw::geom::Extent2D offset(coeff(2*ncoeff*nexp+3*i), coeff(2*ncoeff*nexp+3*i+1));

            // but must be transformed to mm to be applied.
            auto pixelSize = it->second->getPixelSize();
            auto scaling = afw::geom::LinearTransform::makeScaling(pixelSize.getX(), pixelSize.getY());
            offset = scaling(offset);

            afw::cameraGeom::Orientation newOrientation(
                it->second->getOrientation().getFpPosition() + offset,
                it->second->getOrientation().getReferencePoint(),
                it->second->getOrientation().getYaw() + coeff(2*ncoeff*nexp + 3*i + 2) * afw::geom::radians,
                it->second->getOrientation().getPitch(),
                it->second->getOrientation().getRoll());

            afw::cameraGeom::CameraTransformMap::Transforms newTr;

            // Transform from pixels to focal plane has to be recalculated.
            newTr[afw::cameraGeom::FOCAL_PLANE] = boost::make_shared<afw::geom::AffineXYTransform const>
                                                  (newOrientation.makePixelFpTransform(pixelSize));

            // We should not require any other transformations within meas_mosaic.

            // Replace the old Detector with a shifted one.
            it->second = boost::make_shared<afw::cameraGeom::Detector>(
                it->second->getName(),
                it->second->getId(),
                it->second->getType(),
                it->second->getSerial(),
                it->second->getBBox(),
                it->second->getAmpInfoCatalog(),
                newOrientation,
                it->second->getPixelSize(),
                newTr);
		}
	    } else {
		int i = 0;
		for (CcdSet::iterator it = ccdSet.begin(); it != ccdSet.end(); it++, i++) {
            afw::geom::Extent2D offset(coeff(2*ncoeff*nexp+2*i), coeff(2*ncoeff*nexp+2*i+1));
            auto pixelSize = it->second->getPixelSize();
            auto scaling = afw::geom::LinearTransform::makeScaling(pixelSize.getX(), pixelSize.getY());
            offset = scaling(offset);

            afw::cameraGeom::Orientation newOrientation(
                it->second->getOrientation().getFpPosition() + offset,
                it->second->getOrientation().getReferencePoint(),
                it->second->getOrientation().getYaw(),
                it->second->getOrientation().getPitch(),
                it->second->getOrientation().getRoll());

            afw::cameraGeom::CameraTransformMap::Transforms newTr;

            // Transform from pixels to focal plane has to be recalculated.
            newTr[afw::cameraGeom::FOCAL_PLANE] = boost::make_shared<afw::geom::AffineXYTransform const>
                                                  (newOrientation.makePixelFpTransform(pixelSize));

            // We should not require any other transformations within meas_mosaic.

            // Replace the old Detector with a shifted one.
            it->second = boost::make_shared<afw::cameraGeom::Detector>(
                it->second->getName(),
                it->second->getId(),
                it->second->getType(),
                it->second->getSerial(),
                it->second->getBBox(),
                it->second->getAmpInfoCatalog(),
                newOrientation,
                it->second->getPixelSize(),
                newTr);
		}
	    }
	}

	for (int i = 0; i < nMobs; i++) {
	    matchVec[i]->setUV(ccdSet[matchVec[i]->ichip],
			       coeffVec[matchVec[i]->iexp]->x0, coeffVec[matchVec[i]->iexp]->y0);
	    matchVec[i]->setFitVal(coeffVec[matchVec[i]->iexp], p);
	}

	long size0;
	if (solveCcd) {
	    if (allowRotation) {
		size0 = 2*ncoeff*nexp + 3*nchip + 1;
	    } else {
		size0 = 2*ncoeff*nexp + 2*nchip;
	    }
	} else {
	    size0 = 2*ncoeff*nexp;
	}

	for (int i = 0; i < nSobs; i++) {
	    if (sourceVec[i]->jstar != -1) {
		sourceVec[i]->ra  += coeff(size0+2*sourceVec[i]->jstar);
		sourceVec[i]->dec += coeff(size0+2*sourceVec[i]->jstar+1);
		double rac  = coeffVec[sourceVec[i]->iexp]->A;
		double decc = coeffVec[sourceVec[i]->iexp]->D;
		sourceVec[i]->setXiEta(rac, decc);
		sourceVec[i]->setUV(ccdSet[sourceVec[i]->ichip], coeffVec[sourceVec[i]->iexp]->x0, coeffVec[sourceVec[i]->iexp]->y0);
		sourceVec[i]->setFitVal(coeffVec[sourceVec[i]->iexp], p);
	    } else {
		sourceVec[i]->setUV(ccdSet[sourceVec[i]->ichip], coeffVec[sourceVec[i]->iexp]->x0, coeffVec[sourceVec[i]->iexp]->y0);
		sourceVec[i]->setFitVal(coeffVec[sourceVec[i]->iexp], p);
	    }
	}

	if (writeSnapshots) {
	    writeObsVec((snapshotPath / (boost::format("match-iter-%d.fits") % k).str()).native(), matchVec);
	    writeObsVec((snapshotPath / (boost::format("source-iter-%d.fits") % k).str()).native(), sourceVec);
	}

	double chi2 = calcChi2_Star(matchVec, sourceVec, coeffVec, p);
	printf("%dth iteration calcChi2: %e %e\n", (k+1), calcChi2(matchVec, coeffVec, p), chi2);
	printf("%dth iteration matched: %5.3f (arcsec) sources: %5.3f (arcsec)\n",
	       (k+1),
	       sqrt(calcChi2(matchVec, coeffVec, p, true))*3600.0,
	       sqrt(calcChi2(sourceVec, coeffVec, p, true))*3600.0);
	//flagObj2(matchVec, coeffVec, p, 9.0, catRMS);
	//flagObj2(sourceVec, coeffVec, p, 9.0);
	double drm = calcChi2(matchVec, coeffVec, p, true);
	double drs = calcChi2(sourceVec, coeffVec, p, true);
	flagObj2(matchVec, coeffVec, p, 9.0*drm, catRMS);
	flagObj2(sourceVec, coeffVec, p, 9.0*drs);
    }

    std::map<int, Eigen::Matrix2d> cd;
    for (CoeffSet::iterator it = coeffVec.begin(); it != coeffVec.end(); it++) {
	Eigen::Matrix2d c;
	c << it->second->a[0], it->second->a[1], it->second->b[0], it->second->b[1];
	cd.insert(std::map<int, Eigen::Matrix2d>::value_type(it->first, c));
    }
    for (int i = 0; i < nMobs; i++) {
	double CD1_1 = cd[matchVec[i]->iexp](0,0);
	double CD1_2 = cd[matchVec[i]->iexp](0,1);
	double CD2_1 = cd[matchVec[i]->iexp](1,0);
	double CD2_2 = cd[matchVec[i]->iexp](1,1);
	double det = CD1_1 * CD2_2 - CD1_2 * CD2_1;
	matchVec[i]->U = ( matchVec[i]->xi * CD2_2 - matchVec[i]->eta * CD1_2) / det;
	matchVec[i]->V = (-matchVec[i]->xi * CD2_1 + matchVec[i]->eta * CD1_1) / det;
    }
    for (int i = 0; i < nSobs; i++) {
	double CD1_1 = cd[sourceVec[i]->iexp](0,0);
	double CD1_2 = cd[sourceVec[i]->iexp](0,1);
	double CD2_1 = cd[sourceVec[i]->iexp](1,0);
	double CD2_2 = cd[sourceVec[i]->iexp](1,1);
	double det = CD1_1 * CD2_2 - CD1_2 * CD2_1;
	sourceVec[i]->U = ( sourceVec[i]->xi * CD2_2 - sourceVec[i]->eta * CD1_2) / det;
	sourceVec[i]->V = (-sourceVec[i]->xi * CD2_1 + sourceVec[i]->eta * CD1_1) / det;
    }

    for (CoeffSet::iterator it = coeffVec.begin(); it != coeffVec.end(); it++) {
	std::vector<Obs::Ptr> obsVec_sub;
	for (size_t j = 0; j < matchVec.size(); j++) {
	    Obs::Ptr iobs = matchVec[j];
	    if (iobs->iexp == it->first) {
		obsVec_sub.push_back(iobs);
	    }
	}
	for (size_t j = 0; j < sourceVec.size(); j++) {
	    Obs::Ptr iobs = sourceVec[j];
	    if (iobs->iexp == it->first) {
		obsVec_sub.push_back(iobs);
	    }
	}
	Eigen::VectorXd a = solveSIP_P(p, obsVec_sub);
	for (int k = 0; k < p->ncoeff; k++) {
	    it->second->ap[k] = a(k);
	    it->second->bp[k] = a(k+p->ncoeff);
	}
    }

    for (int i = 0; i < nMobs; i++) {
	matchVec[i]->setFitVal2(coeffVec[matchVec[i]->iexp], p);
    }
    for (int i = 0; i < nSobs; i++) {
	sourceVec[i]->setFitVal2(coeffVec[sourceVec[i]->iexp], p);
    }

    return coeffVec;
}

int fact(int n)
{
    if (n == 1 || n == 0) {
	return 1;
    } else {
	return n * fact(n-1);
    }
}

int binomial(int n, int k)
{
    return (fact(n)/(fact(n-k)*fact(k)));
}

Coeff::Ptr
lsst::meas::mosaic::convertCoeff(Coeff::Ptr& coeff, PTR(lsst::afw::cameraGeom::Detector)& ccd)
{
    Poly::Ptr p = Poly::Ptr(new Poly(coeff->p->order));
    Coeff::Ptr newC = Coeff::Ptr(new Coeff(p));

    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double cosYaw = std::cos(getYaw(ccd));
    double sinYaw = std::sin(getYaw(ccd));

    newC->A = coeff->A;
    newC->D = coeff->D;

    // u = cc * u' - ss * v'
    // v = ss * u' + cc * v'
    // u^i * v^j = (cc * u' - ss * v')^i * (ss * u' + cc * v')^j
    //           = \Sigma (i, n) * (cc * u')^n * (-ss * v')^(i-n) *
    //             \Sigma (j, m) * (ss * u')^m * ( cc * v')^(j-m)
    for (int k = 0; k < p->ncoeff; k++) {
	for (int n = 0; n <= xorder[k]; n++) {
	    for (int m = 0; m <= yorder[k]; m++) {
		int i = n + m;
		int j = xorder[k] + yorder[k] - n - m;
		int l = p->getIndex(i, j);
		double C =  binomial(xorder[k], n) *
		            binomial(yorder[k], m) *
		            pow(cosYaw, n) * pow(-sinYaw, xorder[k]-n) *
		            pow(sinYaw, m) * pow( cosYaw, yorder[k]-m);
		newC->a[l] += coeff->a[k] * C;
		newC->b[l] += coeff->b[k] * C;
	    }
	}
    }

    afw::geom::Extent2D off = getCenterInFpPixels(ccd) - getCenterInDetectorPixels(ccd);

    newC->x0 =  (off[0] + coeff->x0) * cosYaw + (off[1] + coeff->y0) * sinYaw;
    newC->y0 = -(off[0] + coeff->x0) * sinYaw + (off[1] + coeff->y0) * cosYaw;

    double a = coeff->a[0];
    double b = coeff->a[1];
    double c = coeff->b[0];
    double d = coeff->b[1];
    double det = a * d - b * c;
    Eigen::Matrix2d cdinv;
    cdinv << d/det, -b/det, -c/det, a/det;
    Eigen::Matrix2d cd2;
    cd2 << newC->a[0], newC->a[1], newC->b[0], newC->b[1];
    Eigen::Matrix2d mat = cdinv * cd2;
    a = mat(0,0);
    b = mat(0,1);
    c = mat(1,0);
    d = mat(1,1);

    Eigen::VectorXd ap = Eigen::VectorXd::Zero(p->ncoeff);
    Eigen::VectorXd bp = Eigen::VectorXd::Zero(p->ncoeff);

    for (int k = 0; k < p->ncoeff; k++) {
	for (int n = 0; n <= xorder[k]; n++) {
	    for (int m = 0; m <= yorder[k]; m++) {
		int i = n + m;
		int j = xorder[k] + yorder[k] - n - m;
		int l = p->getIndex(i, j);
		double C =  binomial(xorder[k], n) *
		            binomial(yorder[k], m) *
		            pow(a, n) * pow(b, xorder[k]-n) *
		            pow(c, m) * pow(d, yorder[k]-m);
		ap(l) += coeff->ap[k] * C;
		bp(l) += coeff->bp[k] * C;
	    }
	}
    }
    ap(0) += a;
    ap(1) += b;
    bp(0) += c;
    bp(1) += d;

    for (int k = 0; k < p->ncoeff; k++) {
	newC->ap[k] =  ap(k) * cosYaw + bp(k) * sinYaw;
	newC->bp[k] = -ap(k) * sinYaw + bp(k) * cosYaw;
    }
    //newC->ap[0] += cosYaw - 1.;
    //newC->ap[1] += sinYaw;
    //newC->bp[0] -= sinYaw;
    //newC->bp[1] += cosYaw - 1.;
    newC->ap[0] -= 1.;
    newC->bp[1] -= 1.;

    return newC;
}

lsst::afw::image::TanWcs::Ptr
lsst::meas::mosaic::wcsFromCoeff(Coeff::Ptr& coeff)
{
    int order = coeff->p->order;

    lsst::afw::geom::PointD crval
	= lsst::afw::geom::Point2D(coeff->A*R2D, coeff->D*R2D);
    lsst::afw::geom::PointD crpix = lsst::afw::geom::Point2D(-coeff->x0, -coeff->y0);

    Eigen::Matrix2d cd;
    cd << coeff->a[0], coeff->a[1], coeff->b[0], coeff->b[1];
    double D = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    //std::cout << cd << std::endl;

    Eigen::MatrixXd sipA = Eigen::MatrixXd::Zero(order+1,order+1);
    Eigen::MatrixXd sipB = Eigen::MatrixXd::Zero(order+1,order+1);
    for (int k = 2; k <= order; k++) {
	for (int i = k; i >= 0; i--) {
	    int j = k - i;
	    int n = k*(k+1)/2 - 1;
	    sipA(i,j) = ( cd(1,1)*coeff->a[n+j] - cd(0,1)*coeff->b[n+j]) / D;
	    sipB(i,j) = (-cd(1,0)*coeff->a[n+j] + cd(0,0)*coeff->b[n+j]) / D;
	}
    }
    //std::cout << sipA << std::endl;
    //std::cout << sipB << std::endl;

    //cd *= R2D;

    Eigen::MatrixXd sipAp = Eigen::MatrixXd::Zero(order+1,order+1);
    Eigen::MatrixXd sipBp = Eigen::MatrixXd::Zero(order+1,order+1);
    for (int k = 1; k <= order; k++) {
	for (int i = k; i >= 0; i--) {
	    int j = k - i;
	    int n = k*(k+1)/2 - 1;
	    sipAp(i,j) = coeff->ap[n+j];
	    sipBp(i,j) = coeff->bp[n+j];
	}
    }
    //std::cout << sipAp << std::endl;
    //std::cout << sipBp << std::endl;

    lsst::afw::image::TanWcs::Ptr wcs = lsst::afw::image::TanWcs::Ptr(new lsst::afw::image::TanWcs(crval, crpix, cd, sipA, sipB, sipAp, sipBp));

    return wcs;
}

// wholesale copied from afw::image::TanWcs.cc
///@brief Decode the SIP headers for a given matrix, if present.
static void decodeSipHeader(CONST_PTR(lsst::daf::base::PropertySet) const& fitsMetadata,
                            std::string const& which,
                            Eigen::MatrixXd *m) {
    std::string header = which + "_ORDER";
    if (!fitsMetadata->exists(header)) return;
    int order = fitsMetadata->getAsInt(header);
    m->resize(order + 1, order + 1);
    boost::format format("%1%_%2%_%3%");
    for (int i = 0; i <= order; ++i) {
        for (int j = 0; j <= order; ++j) {
            header = (format % which % i % j).str();
            if (fitsMetadata->exists(header)) {
                (*m)(i,j) = fitsMetadata->getAsDouble(header);
            }
            else {
                (*m)(i, j) = 0.0;
            }
        }
    }
}

//hsc::meas::mosaic::coeffFromTanWcs(lsst::afw::image::TanWcs::Ptr& tanwcs)
Coeff::Ptr lsst::meas::mosaic::coeffFromTanWcs(lsst::afw::image::Wcs::Ptr& wcs)
{

//    lsst::daf::base::PropertyList::Ptr fitsMetadata = tanwcs->getFitsMetadata();
    lsst::daf::base::PropertyList::Ptr fitsMetadata = wcs->getFitsMetadata();
    int orderA = fitsMetadata->get<int>("A_ORDER");
    int orderB = fitsMetadata->get<int>("B_ORDER");
    int orderAP = fitsMetadata->get<int>("AP_ORDER");
    int orderBP = fitsMetadata->get<int>("BP_ORDER");
    int order = orderA;
    int orderP = orderAP;

    Poly::Ptr p = Poly::Ptr(new Poly( order ));
    Coeff::Ptr coeff = Coeff::Ptr(new Coeff(p));

    double cd00 = fitsMetadata->get<double>("CD1_1");
    double cd11 = fitsMetadata->get<double>("CD2_2");
    double cd01 = fitsMetadata->get<double>("CD1_2");
    double cd10 = fitsMetadata->get<double>("CD2_1");
    double crval0 = fitsMetadata->get<double>("CRVAL1");
    double crval1 = fitsMetadata->get<double>("CRVAL2");
    double crpix0 = fitsMetadata->get<double>("CRPIX1");
    double crpix1 = fitsMetadata->get<double>("CRPIX2");

    Eigen::MatrixXd sipA = Eigen::MatrixXd::Zero(orderA+1,orderA+1);
    Eigen::MatrixXd sipB = Eigen::MatrixXd::Zero(orderB+1,orderB+1);
    Eigen::MatrixXd sipAp = Eigen::MatrixXd::Zero(orderAP+1,orderAP+1);
    Eigen::MatrixXd sipBp = Eigen::MatrixXd::Zero(orderBP+1,orderBP+1);
    decodeSipHeader(fitsMetadata, "A", &sipA);
    decodeSipHeader(fitsMetadata, "B", &sipB);
    decodeSipHeader(fitsMetadata, "AP", &sipAp);
    decodeSipHeader(fitsMetadata, "BP", &sipBp);

    coeff->set_A(crval0*D2R);
    coeff->set_D(crval1*D2R);
    coeff->set_x0(-crpix0);
    coeff->set_y0(-crpix1);

    coeff->set_iexp(0);

    for (int k = 2; k <= order; k++) {
	for (int i = k; i >= 0; i--) {
	    int j = k - i;
	    int n = k*(k+1)/2 - 1;

            std::cout << "sipA(" << i << "," << j << "): " << sipA(i,j) << std::endl;
            std::cout << "sipB(" << i << "," << j << "): " << sipB(i,j) << std::endl;
            coeff->set_a(n+j, cd00*sipA(i,j)+cd01*sipB(i,j));
            coeff->set_b(n+j, cd10*sipA(i,j)+cd11*sipB(i,j));

	}
    }
    //std::cout << sipA << std::endl;
    //std::cout << sipB << std::endl;

    //cd *= R2D;

    for (int k = 1; k <= orderP; k++) {
	for (int i = k; i >= 0; i--) {
	    int j = k - i;
	    int n = k*(k+1)/2 - 1;

	    coeff->set_ap(n+j, sipAp(i,j));
            coeff->set_bp(n+j, sipBp(i,j));
	}
    }
    //std::cout << sipAp << std::endl;
    //std::cout << sipBp << std::endl;

    // setting CD matrix to coeff.a(0,1) & b(0,1)
    coeff->set_a(0, cd00);
    coeff->set_a(1, cd01);
    coeff->set_b(0, cd10);
    coeff->set_b(1, cd11);
    //double D = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    //std::cout << cd << std::endl;

    return coeff;
}

lsst::afw::image::Image<float>::Ptr
lsst::meas::mosaic::getJImg(Coeff::Ptr& coeff,
			   PTR(lsst::afw::cameraGeom::Detector)& ccd)
{
    double scale = coeff->pixelScale();
    double deg2pix = 1. / scale;

    int width  = getWidth(ccd);
    int height = getHeight(ccd);

    lsst::afw::image::Image<float>::Ptr img(new lsst::afw::image::Image<float>(width, height));

    Eigen::VectorXd vals(width);

    int interpLength = 100;

    for (int y = 0; y != height; y++) {

	for (int x = 0; x < width + interpLength; x+= interpLength) {
	    int interval = interpLength;
	    int xend = x + interval - 1;
	    if (xend >= width) {
		xend = width - 1;
		interval = xend - x + 1;
	    }

        afw::geom::Point2D uv = detPxToFpPx(ccd, afw::geom::Point2D(x, y)) +
                                afw::geom::Extent2D(coeff->x0, coeff->y0);
	    double val0 = coeff->detJ(uv.getX(), uv.getY()) * deg2pix * deg2pix;
        uv = detPxToFpPx(ccd, afw::geom::Point2D(xend, y)) + afw::geom::Extent2D(coeff->x0, coeff->y0);
	    double val1 = coeff->detJ(uv.getX(), uv.getY()) * deg2pix * deg2pix;

	    for (int i = 0; i < interval; i++) {
		vals(x+i) = val0 + (val1 - val0) / interval * i;
	    }
	}

	lsst::afw::image::Image<float>::x_iterator begin = img->row_begin(y);
	lsst::afw::image::Image<float>::x_iterator end   = img->row_end(y);

	for (lsst::afw::image::Image<float>::x_iterator ptr = begin; ptr != end; ptr++) {

	    int x = ptr - begin;

	    *ptr = vals(x);
	}
    }

    return img;
}

lsst::afw::image::Image<float>::Ptr
lsst::meas::mosaic::getJImg(lsst::afw::image::Wcs::Ptr& wcs,
			   int width, int height)
{
    double scale = wcs->pixelScale().asDegrees();
    double deg2pix = 1. / scale;

    lsst::afw::image::Image<float>::Ptr img(new lsst::afw::image::Image<float>(width, height));

    Eigen::VectorXd vals(width);

    int interpLength = 100;

    for (int y = 0; y != height; y++) {

	for (int x = 0; x < width + interpLength; x+= interpLength) {
	    int interval = interpLength;
	    int xend = x + interval - 1;
	    if (xend >= width) {
		xend = width - 1;
		interval = xend - x + 1;
	    }

	    double u = x;
	    double v = y;
	    double val0 = wcs->pixArea(lsst::afw::geom::Point2D(u, v)) * deg2pix * deg2pix;
	    u = xend;
	    v = y;
	    double val1 = wcs->pixArea(lsst::afw::geom::Point2D(u, v)) * deg2pix * deg2pix;

	    for (int i = 0; i < interval; i++) {
		vals(x+i) = val0 + (val1 - val0) / interval * i;
	    }
	}

	lsst::afw::image::Image<float>::x_iterator begin = img->row_begin(y);
	lsst::afw::image::Image<float>::x_iterator end   = img->row_end(y);

	for (lsst::afw::image::Image<float>::x_iterator ptr = begin; ptr != end; ptr++) {

	    int x = ptr - begin;

	    *ptr = vals(x);
	}
    }

    return img;
}

lsst::afw::image::Image<float>::Ptr
lsst::meas::mosaic::getJImg(lsst::afw::image::Wcs::Ptr& wcs,
			   PTR(lsst::afw::cameraGeom::Detector)& ccd)
{
    int width  = getWidth(ccd);
    int height = getHeight(ccd);

    return getJImg(wcs, width, height);
}
