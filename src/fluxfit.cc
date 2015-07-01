#include "ndarray.h"
#include "lsst/meas/mosaic/mosaicfit.h"
#include "lsst/meas/mosaic/fluxfit.h"

using namespace lsst::meas::mosaic;

extern Eigen::VectorXd solveMatrix(long size, Eigen::MatrixXd &a_data, Eigen::VectorXd &b_data);
extern int binomial(int n, int k);

FluxFitParams::FluxFitParams(int order_, bool absolute_, bool chebyshev_) :
    order(order_),
    chebyshev(chebyshev_),
    ncoeff((order+1) * (order+2) / 2),
    xorder(new int[ncoeff]),
    yorder(new int[ncoeff]),
    absolute(absolute_),
    coeff(new double[ncoeff]),
    u_max(1.0),
    v_max(1.0),
    x0(0.0),
    y0(0.0)
{
    int k = 0;
    for (int j = 0; j <= order; j++) {
	for (int i = 0; i <= j; i++) {
	    xorder[k] = j - i;
	    yorder[k] = i;
	    coeff[k] = 0.0;
	    k++;
	}
    }
    assert(k == ncoeff);
}

FluxFitParams::FluxFitParams(lsst::daf::base::PropertySet::Ptr& metadata) :
    order(metadata->getAsInt("ORDER")),
    chebyshev(metadata->getAsBool("CHEBYSHEV")),
    ncoeff((order+1) * (order+2) / 2),
    xorder(new int[ncoeff]),
    yorder(new int[ncoeff]),
    absolute(metadata->getAsBool("ABSOLUTE")),
    coeff(new double[ncoeff]),
    u_max(metadata->getAsDouble("U_MAX")),
    v_max(metadata->getAsDouble("V_MAX")),
    x0(metadata->getAsDouble("X0")),
    y0(metadata->getAsDouble("Y0"))
{
    int k = 0;
    for (int j = 0; j <= order; j++) {
	for (int i = 0; i <= j; i++) {
	    xorder[k] = j - i;
	    yorder[k] = i;
	    std::string label = boost::str(boost::format("C_%d_%d") % xorder[k] % yorder[k]);
	    coeff[k] = metadata->getAsDouble(label);
	    k++;
	}
    }
    assert(k == ncoeff);
}


FluxFitParams::~FluxFitParams(void) {
    delete [] xorder;
    delete [] yorder;
    delete [] coeff;
}

FluxFitParams::FluxFitParams(const FluxFitParams &p) {
    if (p.chebyshev == false) {
	this->order = p.order;
	this->absolute = p.absolute;
	this->chebyshev = p.chebyshev;
	this->ncoeff = p.ncoeff;
	this->u_max = p.u_max;
	this->v_max = p.v_max;
	this->x0 = p.x0;
	this->y0 = p.y0;
	this->xorder = new int[this->ncoeff];
	this->yorder = new int[this->ncoeff];
	this->coeff = new double[this->ncoeff];
	for (int i = 0; i < this->ncoeff; i++) {
	    this->xorder[i] = p.xorder[i];
	    this->yorder[i] = p.yorder[i];
	    this->coeff[i] = p.coeff[i];
	}
    } else {
	this->order = p.order;
	this->absolute = p.absolute;
	this->chebyshev = false;
	this->ncoeff = p.ncoeff;
	this->u_max = p.u_max;
	this->v_max = p.v_max;
	this->x0 = p.x0;
	this->y0 = p.y0;

	Chev c(p.order);
	this->xorder = new int[this->ncoeff];
	this->yorder = new int[this->ncoeff];
	this->coeff = new double[this->ncoeff];
	for (int i = 0; i < this->ncoeff; i++) {
	    this->xorder[i] = p.xorder[i];
	    this->yorder[i] = p.yorder[i];
	    this->coeff[i] = 0.0;
	}
	for (int k = 0; k < this->ncoeff; k++) {
	    for (int i = 0; i <= this->xorder[k]; i++) {
		for (int j = 0; j <= this->yorder[k]; j++) {
		    int kk = this->getIndex(i, j);
		    this->coeff[kk] += p.coeff[k] * 
			               c.coeffs[this->xorder[k]][this->xorder[k]-i] *
			               c.coeffs[this->yorder[k]][this->yorder[k]-j];
		}
	    }
	}
    }
}

double FluxFitParams::eval(double u, double v) const {
   double uu = (u + x0) / u_max;
   double vv = (v + y0) / v_max;
   double val = 0.0;

   if (this->chebyshev) {
      for (int k = 0; k < ncoeff; k++) {
	 val += coeff[k] * Tn(xorder[k], uu) * Tn(yorder[k], vv);
      }
   } else {
      for (int k = 0; k < ncoeff; k++) {
	 val += coeff[k] * pow(uu, xorder[k]) * pow(vv, yorder[k]);
      }
   }

   return val;
}


ndarray::Array<double, 1> FluxFitParams::eval(
    ndarray::Array<double const, 1> const& x,
    ndarray::Array<double const, 1> const& y
    ) const
{
    int const num = x.getShape()[0];
    if (y.getShape()[0] != num) {
        throw LSST_EXCEPT(pex::exceptions::LengthError,
                          str(boost::format("Size mismatch: %d vs %d") % x.getShape()[0] % y.getShape()[0]));
    }
    ndarray::Array<double, 1> out = ndarray::allocate(ndarray::makeVector(num));
    for (int i = 0; i < num; ++i) {
        out[i] = eval(x[i], y[i]);
    }
    return out;
}


int FluxFitParams::getIndex(int i, int j) const {
    for (int k = 0; k < this->ncoeff; k++) {
	if (xorder[k] == i &&
	    yorder[k] == j) {
	    return k;
	}
    }

    return -1;
}

Eigen::VectorXd fluxFit_rel(std::vector<Obs::Ptr> &m,
			    int nmatch,
			    std::vector<Obs::Ptr> &s,
			    int nsource,
			    int nexp,
			    int nchip,
			    FfpSet &ffpSet,
			    bool solveCcd)
{
    int nMobs = m.size();
    int nSobs = s.size();
    int nFfp = ffpSet.size();

    int* num = new int[nmatch+nsource];
    for (int i = 0; i < nmatch+nsource; i++) {
	num[i] = 0;
    }
    for (int i = 0; i < nMobs; i++) {
	if (m[i]->good && m[i]->mag != -9999 && m[i]->err != -9999) {
	    num[m[i]->istar] += 1;
	}
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good && s[i]->mag != -9999 && s[i]->err != -9999) {
	    num[nmatch+s[i]->istar] += 1;
	}
    }
    std::vector<int> v_istar;
    for (int i = 0; i < nmatch+nsource; i++) {
	if (num[i] >= 2) {
	    v_istar.push_back(i);
	}
    }
    delete [] num;
    int nstar = v_istar.size();
    std::cout << "nstar: " << nstar << std::endl;

    for (int i = 0; i < nMobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), m[i]->istar);
	if (it != v_istar.end()) {
	    m[i]->jstar = it - v_istar.begin();
	} else {
	    m[i]->jstar = -1;
	}
    }
    for (int i = 0; i < nSobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), nmatch+s[i]->istar);
	if (it != v_istar.end()) {
	    s[i]->jstar = it - v_istar.begin();
	} else {
	    s[i]->jstar = -1;
	}
    }

    int ncoeff_offset = 3;	// Fit from 2nd order only
    // In some cases (small number of visits or small dithering),
    // fitting from 1st order will degenerate and fails.
    // Currently I'm fitting from 2nd order
    FluxFitParams::Ptr p = ffpSet[ffpSet.begin()->first];
    int ncoeff = p->ncoeff - ncoeff_offset;
    int *xorder = &p->xorder[ncoeff_offset];
    int *yorder = &p->yorder[ncoeff_offset];
    double u_max = p->u_max;
    double v_max = p->v_max;

    Eigen::VectorXd pu(ncoeff);
    Eigen::VectorXd pv(ncoeff);

    int ndim;
    if (solveCcd) {
	ndim = nexp + nchip + ncoeff * nFfp + nstar + 2;
    } else {
	ndim = nexp + ncoeff * nFfp + nstar + 1;
    }
    std::cout << "ndim: " << ndim << std::endl;

    Eigen::MatrixXd a_data = Eigen::MatrixXd::Zero(ndim, ndim);
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(ndim);

    double is2 = 1.0;

    if (solveCcd) {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 || m[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], m[i]->u/u_max);
		    pv(k) = Tn(yorder[k], m[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(m[i]->u/u_max, xorder[k]);
		    pv(k) = pow(m[i]->v/v_max, yorder[k]);
		}
	    }
 
	    is2 = 1.0 / pow(m[i]->err, 2);

	    a_data(m[i]->jexp, m[i]->jexp) -= is2;
	    a_data(m[i]->jexp, nexp+m[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(m[i]->jexp, nexp+nchip+ncoeff*m[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(m[i]->jexp, nexp+nchip+ncoeff*nFfp+m[i]->jstar) += is2;

	    a_data(nexp+m[i]->jchip, m[i]->jexp) -= is2;
	    a_data(nexp+m[i]->jchip, nexp+m[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+m[i]->jchip, nexp+nchip+ncoeff*m[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+m[i]->jchip, nexp+nchip+ncoeff*nFfp+m[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+nchip+ncoeff*m[i]->jexp+j, m[i]->jexp) -= pu(j) * pv(j) * is2;
		a_data(nexp+nchip+ncoeff*m[i]->jexp+j, nexp+m[i]->jchip) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+nchip+ncoeff*m[i]->jexp+j, nexp+nchip+ncoeff*m[i]->jexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+nchip+ncoeff*m[i]->jexp+j, nexp+nchip+ncoeff*nFfp+m[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+nchip+ncoeff*nFfp+m[i]->jstar, m[i]->jexp) += is2;
	    a_data(nexp+nchip+ncoeff*nFfp+m[i]->jstar, nexp+m[i]->jchip) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+nchip+ncoeff*nFfp+m[i]->jstar, nexp+nchip+ncoeff*m[i]->jexp+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+nchip+ncoeff*nFfp+m[i]->jstar, nexp+nchip+ncoeff*nFfp+m[i]->jstar) -= is2;

	    b_data(m[i]->jexp) += m[i]->mag * is2;
	    b_data(nexp+m[i]->jchip) += m[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+nchip+ncoeff*m[i]->jexp+k) += m[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+nchip+ncoeff*nFfp+m[i]->jstar) -= m[i]->mag * is2;
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], s[i]->u/u_max);
		    pv(k) = Tn(yorder[k], s[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(s[i]->u/u_max, xorder[k]);
		    pv(k) = pow(s[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / pow(s[i]->err, 2);

	    a_data(s[i]->jexp, s[i]->jexp) -= is2;
	    a_data(s[i]->jexp, nexp+s[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(s[i]->jexp, nexp+nchip+ncoeff*s[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(s[i]->jexp, nexp+nchip+ncoeff*nFfp+s[i]->jstar) += is2;

	    a_data(nexp+s[i]->jchip, s[i]->jexp) -= is2;
	    a_data(nexp+s[i]->jchip, nexp+s[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+s[i]->jchip, nexp+nchip+ncoeff*s[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+s[i]->jchip, nexp+nchip+ncoeff*nFfp+s[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+nchip+ncoeff*s[i]->jexp+j, s[i]->jexp) -= pu(j) * pv(j) * is2;
		a_data(nexp+nchip+ncoeff*s[i]->jexp+j, nexp+s[i]->jchip) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+nchip+ncoeff*s[i]->jexp+j, nexp+nchip+ncoeff*s[i]->jexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+nchip+ncoeff*s[i]->jexp+j, nexp+nchip+ncoeff*nFfp+s[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+nchip+ncoeff*nFfp+s[i]->jstar, s[i]->jexp) += is2;
	    a_data(nexp+nchip+ncoeff*nFfp+s[i]->jstar, nexp+s[i]->jchip) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+nchip+ncoeff*nFfp+s[i]->jstar, nexp+nchip+ncoeff*s[i]->jexp+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+nchip+ncoeff*nFfp+s[i]->jstar, nexp+nchip+ncoeff*nFfp+s[i]->jstar) -= is2;

	    b_data(s[i]->jexp) += s[i]->mag * is2;
	    b_data(nexp+s[i]->jchip) += s[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+nchip+ncoeff*s[i]->jexp+k) += s[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+nchip+ncoeff*nFfp+s[i]->jstar) -= s[i]->mag * is2;
	}

	a_data(0, nexp+nchip+ncoeff*nFfp+nstar) = 1;
	a_data(nexp+nchip+ncoeff*nFfp+nstar, 0) = 1;
	b_data(ndim-2) = 0;

	for (int i = 0; i < nchip; i++) {
	    a_data(nexp+i, nexp+nchip+ncoeff*nFfp+nstar+1) = -1;
	    a_data(nexp+nchip+ncoeff*nFfp+nstar+1, nexp+i) = -1;
	}
	b_data(ndim-1) = 0;
    } else {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 || m[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], m[i]->u/u_max);
		    pv(k) = Tn(yorder[k], m[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(m[i]->u/u_max, xorder[k]);
		    pv(k) = pow(m[i]->v/v_max, yorder[k]);
		}
	    }
 
	    is2 = 1.0 / pow(m[i]->err, 2);

	    a_data(m[i]->jexp, m[i]->jexp) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(m[i]->jexp, nexp+ncoeff*m[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(m[i]->jexp, nexp+ncoeff*nFfp+m[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+ncoeff*m[i]->jexp+j, m[i]->jexp) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+ncoeff*m[i]->jexp+j, nexp+ncoeff*m[i]->jexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+ncoeff*m[i]->jexp+j, nexp+ncoeff*nFfp+m[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+ncoeff*nFfp+m[i]->jstar, m[i]->jexp) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+ncoeff*nFfp+m[i]->jstar, nexp+ncoeff*m[i]->jexp+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+ncoeff*nFfp+m[i]->jstar, nexp+ncoeff*nFfp+m[i]->jstar) -= is2;

	    b_data(m[i]->jexp) += m[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+ncoeff*m[i]->jexp+k) += m[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+ncoeff*nFfp+m[i]->jstar) -= m[i]->mag * is2;
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], s[i]->u/u_max);
		    pv(k) = Tn(yorder[k], s[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(s[i]->u/u_max, xorder[k]);
		    pv(k) = pow(s[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / pow(s[i]->err, 2);

	    a_data(s[i]->jexp, s[i]->jexp) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(s[i]->jexp, nexp+ncoeff*s[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(s[i]->jexp, nexp+ncoeff*nFfp+s[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+ncoeff*s[i]->jexp+j, s[i]->jexp) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+ncoeff*s[i]->jexp+j, nexp+ncoeff*s[i]->jexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+ncoeff*s[i]->jexp+j, nexp+ncoeff*nFfp+s[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+ncoeff*nFfp+s[i]->jstar, s[i]->jexp) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+ncoeff*nFfp+s[i]->jstar, nexp+ncoeff*s[i]->jexp+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+ncoeff*nFfp+s[i]->jstar, nexp+ncoeff*nFfp+s[i]->jstar) -= is2;

	    b_data(s[i]->jexp) += s[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+ncoeff*s[i]->jexp+k) += s[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+ncoeff*nFfp+s[i]->jstar) -= s[i]->mag * is2;
	}

	a_data(0, nexp+ncoeff*nFfp+nstar) = 1;
	a_data(nexp+ncoeff*nFfp+nstar, 0) = 1;
	b_data(ndim-1) = 0;
    }

    Eigen::VectorXd solution = solveMatrix(ndim, a_data, b_data);

    std::vector<double> v;
    std::vector<double> e;
    for (int i = 0; i < nMobs; i++) {
	if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 || m[i]->err == -9999 ||
	    m[i]->mag_cat == -9999) continue;
	if (solveCcd) {
	    v.push_back(m[i]->mag_cat - solution(nexp+nchip+ncoeff*nFfp+m[i]->jstar));
	} else {
	    v.push_back(m[i]->mag_cat - solution(nexp+ncoeff*nFfp+m[i]->jstar));
	}
	e.push_back(m[i]->err_cat);
    }

    double S = 0.;
    double Sx = 0.;
    double Sxx = 0.;
    for (unsigned int i = 0; i < v.size(); i++) {
	S += 1./(e[i]*e[i]);
	Sx += v[i]/(e[i]*e[i]);
	Sxx += v[i]*v[i]/(e[i]*e[i]);
    }
    double avg = Sx / S;
    double std = sqrt((Sxx-Sx*Sx/S)/S);
    std::cout << avg << " " << std << std::endl;

    for (int k = 0; k < 2; k++) {
	S = Sx = Sxx = 0.;
	for (unsigned int i = 0; i < v.size(); i++) {
	    if (fabs(v[i]-avg)/e[i] < 3.0) {
		S += 1./(e[i]*e[i]);
		Sx += v[i]/(e[i]*e[i]);
		Sxx += v[i]*v[i]/(e[i]*e[i]);
	    }
	}
	avg = Sx / S;
	std = sqrt((Sxx-Sx*Sx/S)/S);
	std::cout << avg << " " << std << std::endl;
    }

    double dmag = avg;

    for (int i = 0; i < nexp; i++) {
	solution(i) += dmag;
    }
    if (solveCcd) {
	for (int i = 0; i < nstar; i++) {
	    solution(nexp+nchip+ncoeff*nFfp+i) += dmag;
	}
    } else {
	for (int i = 0; i < nstar; i++) {
	    solution(nexp+ncoeff*nFfp+i) += dmag;
	}
    }

    if (solveCcd) {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999) continue;
	    m[i]->mag0 = solution(nexp+nchip+ncoeff*nFfp+m[i]->jstar);
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	    s[i]->mag0 = solution(nexp+nchip+ncoeff*nFfp+s[i]->jstar);
	}
    } else {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999) continue;
	    m[i]->mag0 = solution(nexp+ncoeff*nFfp+m[i]->jstar);
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	    s[i]->mag0 = solution(nexp+ncoeff*nFfp+s[i]->jstar);
	}
    }

    if (solveCcd) {
	int j = 0;
	for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++, j++) {
	    for (int i = 0; i < ncoeff; i++) {
		ffpSet[it->first]->coeff[ncoeff_offset+i] = solution(nexp+nchip+ncoeff*j+i);
	    }
	}
    } else {
	int j = 0;
	for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++, j++) {
	    for (int i = 0; i < ncoeff; i++) {
		ffpSet[it->first]->coeff[ncoeff_offset+i] = solution(nexp+ncoeff*j+i);
	    }
	}
    }

    return solution;
}

Eigen::VectorXd fluxFit_rel1(std::vector<Obs::Ptr> &m,
			     int nmatch,
			     std::vector<Obs::Ptr> &s,
			     int nsource,
			     int nexp,
			     int nchip,
			     FfpSet &ffpSet,
			     bool solveCcd)
{
    int nMobs = m.size();
    int nSobs = s.size();

    int* num = new int[nmatch+nsource];
    for (int i = 0; i < nmatch+nsource; i++) {
	num[i] = 0;
    }
    for (int i = 0; i < nMobs; i++) {
	if (m[i]->good && m[i]->mag != -9999 && m[i]->err != -9999) {
	    num[m[i]->istar] += 1;
	}
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good && s[i]->mag != -9999 && s[i]->err != -9999) {
	    num[nmatch+s[i]->istar] += 1;
	}
    }
    std::vector<int> v_istar;
    for (int i = 0; i < nmatch+nsource; i++) {
	if (num[i] >= 2) {
	    v_istar.push_back(i);
	}
    }
    delete [] num;
    int nstar = v_istar.size();
    std::cout << "nstar: " << nstar << std::endl;

    for (int i = 0; i < nMobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), m[i]->istar);
	if (it != v_istar.end()) {
	    m[i]->jstar = it - v_istar.begin();
	} else {
	    m[i]->jstar = -1;
	}
    }
    for (int i = 0; i < nSobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), nmatch+s[i]->istar);
	if (it != v_istar.end()) {
	    s[i]->jstar = it - v_istar.begin();
	} else {
	    s[i]->jstar = -1;
	}
    }

    int ncoeff_offset = 3;	// Fit from 2nd order only
    // In some cases (small number of visits or small dithering),
    // fitting from 1st order will degenerate and fails.
    // Currently I'm fitting from 2nd order
    FluxFitParams::Ptr p = ffpSet[ffpSet.begin()->first];
    int ncoeff = p->ncoeff - ncoeff_offset;
    int *xorder = &p->xorder[ncoeff_offset];
    int *yorder = &p->yorder[ncoeff_offset];
    double u_max = p->u_max;
    double v_max = p->v_max;

    Eigen::VectorXd pu(ncoeff);
    Eigen::VectorXd pv(ncoeff);

    int ndim;
    if (solveCcd) {
	ndim = nexp + nchip + ncoeff + nstar + 2;
    } else {
	ndim = nexp + ncoeff + nstar + 1;
    }
    std::cout << "ndim: " << ndim << std::endl;

    Eigen::MatrixXd a_data = Eigen::MatrixXd::Zero(ndim, ndim);
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(ndim);

    double is2 = 1.0;

    if (solveCcd) {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 || m[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], m[i]->u/u_max);
		    pv(k) = Tn(yorder[k], m[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(m[i]->u/u_max, xorder[k]);
		    pv(k) = pow(m[i]->v/v_max, yorder[k]);
		}
	    }
 
	    is2 = 1.0 / pow(m[i]->err, 2);

	    a_data(m[i]->jexp, m[i]->jexp) -= is2;
	    a_data(m[i]->jexp, nexp+m[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(m[i]->jexp, nexp+nchip+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(m[i]->jexp, nexp+nchip+ncoeff+m[i]->jstar) += is2;

	    a_data(nexp+m[i]->jchip, m[i]->jexp) -= is2;
	    a_data(nexp+m[i]->jchip, nexp+m[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+m[i]->jchip, nexp+nchip+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+m[i]->jchip, nexp+nchip+ncoeff+m[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+nchip+j, m[i]->jexp) -= pu(j) * pv(j) * is2;
		a_data(nexp+nchip+j, nexp+m[i]->jchip) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+nchip+j, nexp+nchip+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+nchip+j, nexp+nchip+ncoeff+m[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+nchip+ncoeff+m[i]->jstar, m[i]->jexp) += is2;
	    a_data(nexp+nchip+ncoeff+m[i]->jstar, nexp+m[i]->jchip) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+nchip+ncoeff+m[i]->jstar, nexp+nchip+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+nchip+ncoeff+m[i]->jstar, nexp+nchip+ncoeff+m[i]->jstar) -= is2;

	    b_data(m[i]->jexp) += m[i]->mag * is2;
	    b_data(nexp+m[i]->jchip) += m[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+nchip+k) += m[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+nchip+ncoeff+m[i]->jstar) -= m[i]->mag * is2;
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], s[i]->u/u_max);
		    pv(k) = Tn(yorder[k], s[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(s[i]->u/u_max, xorder[k]);
		    pv(k) = pow(s[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / pow(s[i]->err, 2);

	    a_data(s[i]->jexp, s[i]->jexp) -= is2;
	    a_data(s[i]->jexp, nexp+s[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(s[i]->jexp, nexp+nchip+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(s[i]->jexp, nexp+nchip+ncoeff+s[i]->jstar) += is2;

	    a_data(nexp+s[i]->jchip, s[i]->jexp) -= is2;
	    a_data(nexp+s[i]->jchip, nexp+s[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+s[i]->jchip, nexp+nchip+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+s[i]->jchip, nexp+nchip+ncoeff+s[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+nchip+j, s[i]->jexp) -= pu(j) * pv(j) * is2;
		a_data(nexp+nchip+j, nexp+s[i]->jchip) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+nchip+j, nexp+nchip+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+nchip+j, nexp+nchip+ncoeff+s[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+nchip+ncoeff+s[i]->jstar, s[i]->jexp) += is2;
	    a_data(nexp+nchip+ncoeff+s[i]->jstar, nexp+s[i]->jchip) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+nchip+ncoeff+s[i]->jstar, nexp+nchip+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+nchip+ncoeff+s[i]->jstar, nexp+nchip+ncoeff+s[i]->jstar) -= is2;

	    b_data(s[i]->jexp) += s[i]->mag * is2;
	    b_data(nexp+s[i]->jchip) += s[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+nchip+k) += s[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+nchip+ncoeff+s[i]->jstar) -= s[i]->mag * is2;
	}

	a_data(0, nexp+nchip+ncoeff+nstar) = 1;
	a_data(nexp+nchip+ncoeff+nstar, 0) = 1;
	b_data(ndim-2) = 0;

	for (int i = 0; i < nchip; i++) {
	    a_data(nexp+i, nexp+nchip+ncoeff+nstar+1) = -1;
	    a_data(nexp+nchip+ncoeff+nstar+1, nexp+i) = -1;
	}
	b_data(ndim-1) = 0;
    } else {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 || m[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], m[i]->u/u_max);
		    pv(k) = Tn(yorder[k], m[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(m[i]->u/u_max, xorder[k]);
		    pv(k) = pow(m[i]->v/v_max, yorder[k]);
		}
	    }
 
	    is2 = 1.0 / pow(m[i]->err, 2);

	    a_data(m[i]->jexp, m[i]->jexp) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(m[i]->jexp, nexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(m[i]->jexp, nexp+ncoeff+m[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+j, m[i]->jexp) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+j, nexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+j, nexp+ncoeff+m[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+ncoeff+m[i]->jstar, m[i]->jexp) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+ncoeff+m[i]->jstar, nexp+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+ncoeff+m[i]->jstar, nexp+ncoeff+m[i]->jstar) -= is2;

	    b_data(m[i]->jexp) += m[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+k) += m[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+ncoeff+m[i]->jstar) -= m[i]->mag * is2;
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], s[i]->u/u_max);
		    pv(k) = Tn(yorder[k], s[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(s[i]->u/u_max, xorder[k]);
		    pv(k) = pow(s[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / pow(s[i]->err, 2);

	    a_data(s[i]->jexp, s[i]->jexp) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(s[i]->jexp, nexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(s[i]->jexp, nexp+ncoeff+s[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+j, s[i]->jexp) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+j, nexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+j, nexp+ncoeff+s[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+ncoeff+s[i]->jstar, s[i]->jexp) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+ncoeff+s[i]->jstar, nexp+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+ncoeff+s[i]->jstar, nexp+ncoeff+s[i]->jstar) -= is2;

	    b_data(s[i]->jexp) += s[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+k) += s[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+ncoeff+s[i]->jstar) -= s[i]->mag * is2;
	}

	a_data(0, nexp+ncoeff+nstar) = 1;
	a_data(nexp+ncoeff+nstar, 0) = 1;
	b_data(ndim-1) = 0;
    }

    Eigen::VectorXd solution = solveMatrix(ndim, a_data, b_data);

    std::vector<double> v;
    std::vector<double> e;
    for (int i = 0; i < nMobs; i++) {
	if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 || m[i]->err == -9999 ||
	    m[i]->mag_cat == -9999) continue;
	if (solveCcd) {
	    v.push_back(m[i]->mag_cat - solution(nexp+nchip+ncoeff+m[i]->jstar));
	} else {
	    v.push_back(m[i]->mag_cat - solution(nexp+ncoeff+m[i]->jstar));
	}
	e.push_back(m[i]->err_cat);
    }

    double S = 0.;
    double Sx = 0.;
    double Sxx = 0.;
    for (unsigned int i = 0; i < v.size(); i++) {
	S += 1./(e[i]*e[i]);
	Sx += v[i]/(e[i]*e[i]);
	Sxx += v[i]*v[i]/(e[i]*e[i]);
    }
    double avg = Sx / S;
    double std = sqrt((Sxx-Sx*Sx/S)/S);
    std::cout << avg << " " << std << std::endl;

    for (int k = 0; k < 2; k++) {
	S = Sx = Sxx = 0.;
	for (unsigned int i = 0; i < v.size(); i++) {
	    if (fabs(v[i]-avg)/e[i] < 3.0) {
		S += 1./(e[i]*e[i]);
		Sx += v[i]/(e[i]*e[i]);
		Sxx += v[i]*v[i]/(e[i]*e[i]);
	    }
	}
	avg = Sx / S;
	std = sqrt((Sxx-Sx*Sx/S)/S);
	std::cout << avg << " " << std << std::endl;
    }

    double dmag = avg;

    for (int i = 0; i < nexp; i++) {
	solution(i) += dmag;
    }
    if (solveCcd) {
	for (int i = 0; i < nstar; i++) {
	    solution(nexp+nchip+ncoeff+i) += dmag;
	}
    } else {
	for (int i = 0; i < nstar; i++) {
	    solution(nexp+ncoeff+i) += dmag;
	}
    }

    if (solveCcd) {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999) continue;
	    m[i]->mag0 = solution(nexp+nchip+ncoeff+m[i]->jstar);
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	    s[i]->mag0 = solution(nexp+nchip+ncoeff+s[i]->jstar);
	}
    } else {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999) continue;
	    m[i]->mag0 = solution(nexp+ncoeff+m[i]->jstar);
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	    s[i]->mag0 = solution(nexp+ncoeff+s[i]->jstar);
	}
    }

    if (solveCcd) {
	for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++) {
	    for (int i = 0; i < ncoeff; i++) {
		ffpSet[it->first]->coeff[ncoeff_offset+i] = solution(nexp+nchip+i);
	    }
	}
    } else {
	for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++) {
	    for (int i = 0; i < ncoeff; i++) {
		ffpSet[it->first]->coeff[ncoeff_offset+i] = solution(nexp+i);
	    }
	}
    }

    return solution;
}

Eigen::VectorXd fluxFit_abs(std::vector<Obs::Ptr> &m,
			    int nmatch,
			    std::vector<Obs::Ptr> &s,
			    int nsource,
			    int nexp,
			    int nchip,
			    FfpSet &ffpSet,
			    bool solveCcd)
{
    int nMobs = m.size();
    int nSobs = s.size();
    int nFfp = ffpSet.size();

    int* num = new int[nsource];
    for (int i = 0; i < nsource; i++) {
	num[i] = 0;
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good && s[i]->mag != -9999 && s[i]->err != -9999) {
	    num[s[i]->istar] += 1;
	}
    }
    std::vector<int> v_istar;
    for (int i = 0; i < nsource; i++) {
	if (num[i] >= 2) {
	    v_istar.push_back(i);
	}
    }
    delete [] num;
    int nstar = v_istar.size();
    std::cout << "nstar: " << nstar << std::endl;

    for (int i = 0; i < nSobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), s[i]->istar);
	if (it != v_istar.end()) {
	    s[i]->jstar = it - v_istar.begin();
	} else {
	    s[i]->jstar = -1;
	}
    }

    int ncoeff_offset = 1;	// Fit from 1st order only
    FluxFitParams::Ptr p = ffpSet[ffpSet.begin()->first];
    int ncoeff = p->ncoeff - ncoeff_offset;
    int *xorder = &p->xorder[ncoeff_offset];
    int *yorder = &p->yorder[ncoeff_offset];
    double u_max = p->u_max;
    double v_max = p->v_max;

    Eigen::VectorXd pu(ncoeff);
    Eigen::VectorXd pv(ncoeff);

    int ndim;
    if (solveCcd) {
	ndim = nexp + nchip + ncoeff * nFfp + nstar + 1;
    } else {
	ndim = nexp + ncoeff * nFfp + nstar;
    }
    std::cout << "ndim: " << ndim << std::endl;

    Eigen::MatrixXd a_data = Eigen::MatrixXd::Zero(ndim, ndim);
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(ndim);

    double is2 = 1.0;

    if (solveCcd) {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 ||
		m[i]->err == -9999 || m[i]->mag_cat == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], m[i]->u/u_max);
		    pv(k) = Tn(yorder[k], m[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(m[i]->u/u_max, xorder[k]);
		    pv(k) = pow(m[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / (pow(m[i]->err, 2) + pow(m[i]->err_cat, 2));

	    a_data(m[i]->jexp, m[i]->jexp) -= is2;
	    a_data(m[i]->jexp, nexp+m[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(m[i]->jexp, nexp+nchip+ncoeff*m[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }

	    a_data(nexp+m[i]->jchip, m[i]->jexp) -= is2;
	    a_data(nexp+m[i]->jchip, nexp+m[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+m[i]->jchip, nexp+nchip+ncoeff*m[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+nchip+ncoeff*m[i]->jexp+j, m[i]->jexp) -= pu(j) * pv(j) * is2;
		a_data(nexp+nchip+ncoeff*m[i]->jexp+j, nexp+m[i]->jchip) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+nchip+ncoeff*m[i]->jexp+j, nexp+nchip+ncoeff*m[i]->jexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
	    }

	    b_data(m[i]->jexp) += (m[i]->mag - m[i]->mag_cat) * is2;
	    b_data(nexp+m[i]->jchip) += (m[i]->mag - m[i]->mag_cat) * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+nchip+ncoeff*m[i]->jexp+k) += (m[i]->mag - m[i]->mag_cat) * pu(k) * pv(k) * is2;
	    }
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], s[i]->u/u_max);
		    pv(k) = Tn(yorder[k], s[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(s[i]->u/u_max, xorder[k]);
		    pv(k) = pow(s[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / pow(s[i]->err, 2);

	    a_data(s[i]->jexp, s[i]->jexp) -= is2;
	    a_data(s[i]->jexp, nexp+s[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(s[i]->jexp, nexp+nchip+ncoeff*s[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(s[i]->jexp, nexp+nchip+ncoeff*nFfp+s[i]->jstar) += is2;

	    a_data(nexp+s[i]->jchip, s[i]->jexp) -= is2;
	    a_data(nexp+s[i]->jchip, nexp+s[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+s[i]->jchip, nexp+nchip+ncoeff*s[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+s[i]->jchip, nexp+nchip+ncoeff*nFfp+s[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+nchip+ncoeff*s[i]->jexp+j, s[i]->jexp) -= pu(j) * pv(j) * is2;
		a_data(nexp+nchip+ncoeff*s[i]->jexp+j, nexp+s[i]->jchip) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+nchip+ncoeff*s[i]->jexp+j, nexp+nchip+ncoeff*s[i]->jexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+nchip+ncoeff*s[i]->jexp+j, nexp+nchip+ncoeff*nFfp+s[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+nchip+ncoeff*nFfp+s[i]->jstar, s[i]->jexp) += is2;
	    a_data(nexp+nchip+ncoeff*nFfp+s[i]->jstar, nexp+s[i]->jchip) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+nchip+ncoeff*nFfp+s[i]->jstar, nexp+nchip+ncoeff*s[i]->jexp+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+nchip+ncoeff*nFfp+s[i]->jstar, nexp+nchip+ncoeff*nFfp+s[i]->jstar) -= is2;

	    b_data(s[i]->jexp) += s[i]->mag * is2;
	    b_data(nexp+s[i]->jchip) += s[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+nchip+ncoeff*s[i]->jexp+k) += s[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+nchip+ncoeff*nFfp+s[i]->jstar) -= s[i]->mag * is2;
	}

	for (int i = 0; i < nchip; i++) {
	    a_data(nexp+i, nexp+nchip+ncoeff*nFfp+nstar) = -1;
	    a_data(nexp+nchip+ncoeff*nFfp+nstar, nexp+i) = -1;
	}

	b_data(ndim-1) = 0;
    } else {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 ||
		m[i]->err == -9999 || m[i]->mag_cat == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], m[i]->u/u_max);
		    pv(k) = Tn(yorder[k], m[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(m[i]->u/u_max, xorder[k]);
		    pv(k) = pow(m[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / (pow(m[i]->err, 2) + pow(m[i]->err_cat, 2));

	    a_data(m[i]->jexp, m[i]->jexp) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(m[i]->jexp, nexp+ncoeff*m[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+ncoeff*m[i]->jexp+j, m[i]->jexp) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+ncoeff*m[i]->jexp+j, nexp+ncoeff*m[i]->jexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
	    }

	    b_data(m[i]->jexp) += (m[i]->mag - m[i]->mag_cat) * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+ncoeff*m[i]->jexp+k) += (m[i]->mag - m[i]->mag_cat) * pu(k) * pv(k) * is2;
	    }
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], s[i]->u/u_max);
		    pv(k) = Tn(yorder[k], s[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(s[i]->u/u_max, xorder[k]);
		    pv(k) = pow(s[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / pow(s[i]->err, 2);

	    a_data(s[i]->jexp, s[i]->jexp) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(s[i]->jexp, nexp+ncoeff*s[i]->jexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(s[i]->jexp, nexp+ncoeff*nFfp+s[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+ncoeff*s[i]->jexp+j, s[i]->jexp) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+ncoeff*s[i]->jexp+j, nexp+ncoeff*s[i]->jexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+ncoeff*s[i]->jexp+j, nexp+ncoeff*nFfp+s[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+ncoeff*nFfp+s[i]->jstar, s[i]->jexp) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+ncoeff*nFfp+s[i]->jstar, nexp+ncoeff*s[i]->jexp+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+ncoeff*nFfp+s[i]->jstar, nexp+ncoeff*nFfp+s[i]->jstar) -= is2;

	    b_data(s[i]->jexp) += s[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+ncoeff*s[i]->jexp+k) += s[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+ncoeff*nFfp+s[i]->jstar) -= s[i]->mag * is2;
	}
    }

    Eigen::VectorXd solution = solveMatrix(ndim, a_data, b_data);

    if (solveCcd) {
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	    s[i]->mag0 = solution(nexp+nchip+ncoeff*nFfp+s[i]->jstar);
	}
    } else {
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	    s[i]->mag0 = solution(nexp+ncoeff*nFfp+s[i]->jstar);
	}
    }

    if (solveCcd) {
	int j = 0;
	for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++, j++) {
	    for (int i = 0; i < ncoeff; i++) {
		ffpSet[it->first]->coeff[ncoeff_offset+i] = solution(nexp+nchip+ncoeff*j+i);
	    }
	}
    } else {
	int j = 0;
	for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++, j++) {
	    for (int i = 0; i < ncoeff; i++) {
		ffpSet[it->first]->coeff[ncoeff_offset+i] = solution(nexp+ncoeff*j+i);
	    }
	}
    }

    return solution;
}

Eigen::VectorXd fluxFit_abs1(std::vector<Obs::Ptr> &m,
			     int nmatch,
			     std::vector<Obs::Ptr> &s,
			     int nsource,
			     int nexp,
			     int nchip,
			     //FluxFitParams::Ptr p,
			     FfpSet &ffpSet,
			     bool solveCcd)
{
    int nMobs = m.size();
    int nSobs = s.size();

    int* num = new int[nsource];
    for (int i = 0; i < nsource; i++) {
	num[i] = 0;
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good && s[i]->mag != -9999 && s[i]->err != -9999) {
	    num[s[i]->istar] += 1;
	}
    }
    std::vector<int> v_istar;
    for (int i = 0; i < nsource; i++) {
	if (num[i] >= 2) {
	    v_istar.push_back(i);
	}
    }
    delete [] num;
    int nstar = v_istar.size();
    std::cout << "nstar: " << nstar << std::endl;

    for (int i = 0; i < nSobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), s[i]->istar);
	if (it != v_istar.end()) {
	    s[i]->jstar = it - v_istar.begin();
	} else {
	    s[i]->jstar = -1;
	}
    }

    int ncoeff_offset = 1;	// Fit from 1st order only
    FluxFitParams::Ptr p = ffpSet[ffpSet.begin()->first];
    int ncoeff = p->ncoeff - ncoeff_offset;
    int *xorder = &p->xorder[ncoeff_offset];
    int *yorder = &p->yorder[ncoeff_offset];
    double u_max = p->u_max;
    double v_max = p->v_max;

    Eigen::VectorXd pu(ncoeff);
    Eigen::VectorXd pv(ncoeff);

    int ndim;
    if (solveCcd) {
	ndim = nexp + nchip + ncoeff + nstar + 1;
    } else {
	ndim = nexp + ncoeff + nstar;
    }
    std::cout << "ndim: " << ndim << std::endl;

    Eigen::MatrixXd a_data = Eigen::MatrixXd::Zero(ndim, ndim);
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(ndim);

    double is2 = 1.0;

    if (solveCcd) {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 ||
		m[i]->err == -9999 || m[i]->mag_cat == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], m[i]->u/u_max);
		    pv(k) = Tn(yorder[k], m[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(m[i]->u/u_max, xorder[k]);
		    pv(k) = pow(m[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / (pow(m[i]->err, 2) + pow(m[i]->err_cat, 2));

	    a_data(m[i]->jexp, m[i]->jexp) -= is2;
	    a_data(m[i]->jexp, nexp+m[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(m[i]->jexp, nexp+nchip+k) -= pu(k) * pv(k) * is2;
	    }

	    a_data(nexp+m[i]->jchip, m[i]->jexp) -= is2;
	    a_data(nexp+m[i]->jchip, nexp+m[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+m[i]->jchip, nexp+nchip+k) -= pu(k) * pv(k) * is2;
	    }

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+nchip+j, m[i]->jexp) -= pu(j) * pv(j) * is2;
		a_data(nexp+nchip+j, nexp+m[i]->jchip) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+nchip+j, nexp+nchip+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
	    }

	    b_data(m[i]->jexp) += (m[i]->mag - m[i]->mag_cat) * is2;
	    b_data(nexp+m[i]->jchip) += (m[i]->mag - m[i]->mag_cat) * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+nchip+k) += (m[i]->mag - m[i]->mag_cat) * pu(k) * pv(k) * is2;
	    }
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], s[i]->u/u_max);
		    pv(k) = Tn(yorder[k], s[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(s[i]->u/u_max, xorder[k]);
		    pv(k) = pow(s[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / pow(s[i]->err, 2);

	    a_data(s[i]->jexp, s[i]->jexp) -= is2;
	    a_data(s[i]->jexp, nexp+s[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(s[i]->jexp, nexp+nchip+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(s[i]->jexp, nexp+nchip+ncoeff+s[i]->jstar) += is2;

	    a_data(nexp+s[i]->jchip, s[i]->jexp) -= is2;
	    a_data(nexp+s[i]->jchip, nexp+s[i]->jchip) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+s[i]->jchip, nexp+nchip+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+s[i]->jchip, nexp+nchip+ncoeff+s[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+nchip+j, s[i]->jexp) -= pu(j) * pv(j) * is2;
		a_data(nexp+nchip+j, nexp+s[i]->jchip) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+nchip+j, nexp+nchip+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+nchip+j, nexp+nchip+ncoeff+s[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+nchip+ncoeff+s[i]->jstar, s[i]->jexp) += is2;
	    a_data(nexp+nchip+ncoeff+s[i]->jstar, nexp+s[i]->jchip) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+nchip+ncoeff+s[i]->jstar, nexp+nchip+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+nchip+ncoeff+s[i]->jstar, nexp+nchip+ncoeff+s[i]->jstar) -= is2;

	    b_data(s[i]->jexp) += s[i]->mag * is2;
	    b_data(nexp+s[i]->jchip) += s[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+nchip+k) += s[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+nchip+ncoeff+s[i]->jstar) -= s[i]->mag * is2;
	}

	for (int i = 0; i < nchip; i++) {
	    a_data(nexp+i, nexp+nchip+ncoeff+nstar) = -1;
	    a_data(nexp+nchip+ncoeff+nstar, nexp+i) = -1;
	}

	b_data(ndim-1) = 0;
    } else {
	for (int i = 0; i < nMobs; i++) {
	    if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 ||
		m[i]->err == -9999 || m[i]->mag_cat == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], m[i]->u/u_max);
		    pv(k) = Tn(yorder[k], m[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(m[i]->u/u_max, xorder[k]);
		    pv(k) = pow(m[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / (pow(m[i]->err, 2) + pow(m[i]->err_cat, 2));

	    a_data(m[i]->jexp, m[i]->jexp) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(m[i]->jexp, nexp+k) -= pu(k) * pv(k) * is2;
	    }

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+j, m[i]->jexp) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+j, nexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
	    }

	    b_data(m[i]->jexp) += (m[i]->mag - m[i]->mag_cat) * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+k) += (m[i]->mag - m[i]->mag_cat) * pu(k) * pv(k) * is2;
	    }
	}
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;

	    if (p->chebyshev) {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = Tn(xorder[k], s[i]->u/u_max);
		    pv(k) = Tn(yorder[k], s[i]->v/v_max);
		}
	    } else {
		for (int k = 0; k < ncoeff; k++) {
		    pu(k) = pow(s[i]->u/u_max, xorder[k]);
		    pv(k) = pow(s[i]->v/v_max, yorder[k]);
		}
	    }

	    is2 = 1.0 / pow(s[i]->err, 2);

	    a_data(s[i]->jexp, s[i]->jexp) -= is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(s[i]->jexp, nexp+k) -= pu(k) * pv(k) * is2;
	    }
	    a_data(s[i]->jexp, nexp+ncoeff+s[i]->jstar) += is2;

	    for (int j = 0; j < ncoeff; j++) {
		a_data(nexp+j, s[i]->jexp) -= pu(j) * pv(j) * is2;
		for (int k = 0; k < ncoeff; k++) {
		    a_data(nexp+j, nexp+k) -= pu(j) * pv(j) * pu(k) * pv(k) * is2;
		}
		a_data(nexp+j, nexp+ncoeff+s[i]->jstar) += pu(j) * pv(j) * is2;
	    }

	    a_data(nexp+ncoeff+s[i]->jstar, s[i]->jexp) += is2;
	    for (int k = 0; k < ncoeff; k++) {
		a_data(nexp+ncoeff+s[i]->jstar, nexp+k) += pu(k) * pv(k) * is2;
	    }
	    a_data(nexp+ncoeff+s[i]->jstar, nexp+ncoeff+s[i]->jstar) -= is2;

	    b_data(s[i]->jexp) += s[i]->mag * is2;
	    for (int k = 0; k < ncoeff; k++) {
		b_data(nexp+k) += s[i]->mag * pu(k) * pv(k) * is2;
	    }
	    b_data(nexp+ncoeff+s[i]->jstar) -= s[i]->mag * is2;
	}
    }

    Eigen::VectorXd solution = solveMatrix(ndim, a_data, b_data);

    if (solveCcd) {
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	    s[i]->mag0 = solution(nexp+nchip+ncoeff+s[i]->jstar);
	}
    } else {
	for (int i = 0; i < nSobs; i++) {
	    if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	    s[i]->mag0 = solution(nexp+ncoeff+s[i]->jstar);
	}
    }

    if (solveCcd) {
	for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++) {
	    for (int i = 0; i < ncoeff; i++) {
		ffpSet[it->first]->coeff[ncoeff_offset+i] = solution(nexp+nchip+i);
	    }
	}
    } else {
	for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++) {
	    for (int i = 0; i < ncoeff; i++) {
		ffpSet[it->first]->coeff[ncoeff_offset+i] = solution(nexp+i);
	    }
	}
    }

    return solution;
}

double calcChi2_rel(std::vector<Obs::Ptr> &m, 
		    std::vector<Obs::Ptr> &s,
		    std::map<int, float>& fexp,
		    std::map<int, float>& fchip,
		    //FluxFitParams::Ptr p,
		    FfpSet &ffpSet,
		    bool mag=false)
{
    int nMobs = m.size();
    int nSobs = s.size();

    double chi2 = 0.0;
    double mag2 = 0.0;
    int num = 0;
    for (int i = 0; i < nMobs; i++) {
	if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 || m[i]->err == -9999) continue;
	double val = m[i]->mag - 2.5 * log10(fexp[m[i]->iexp] * fchip[m[i]->ichip]);
	val += ffpSet[m[i]->iexp]->eval(m[i]->u, m[i]->v);
	chi2 += pow((val - m[i]->mag0)/m[i]->err, 2.0);
	mag2 += pow((val - m[i]->mag0), 2.0);
	num++;
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;
	double val = s[i]->mag - 2.5 * log10(fexp[s[i]->iexp] * fchip[s[i]->ichip]);
	val += ffpSet[s[i]->iexp]->eval(s[i]->u, s[i]->v);
	chi2 += pow((val - s[i]->mag0)/s[i]->err, 2.0);
	mag2 += pow((val - s[i]->mag0), 2.0);
	num++;
    }

    if (mag)
	return mag2 / num;
    else
	return chi2 / num;
}

double calcChi2_abs(std::vector<Obs::Ptr> &m, 
		    std::vector<Obs::Ptr> &s,
		    std::map<int, float>& fexp,
		    std::map<int, float>& fchip,
		    //FluxFitParams::Ptr p,
		    FfpSet &ffpSet,
		    bool mag=false)
{
    int nMobs = m.size();
    int nSobs = s.size();

    double chi2 = 0.0;
    double mag2 = 0.0;
    int num = 0;
    for (int i = 0; i < nMobs; i++) {
	if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 ||
	    m[i]->err == -9999 || m[i]->mag_cat == -9999) continue;
	double val = m[i]->mag - 2.5 * log10(fexp[m[i]->iexp] * fchip[m[i]->ichip]);
	val += ffpSet[m[i]->iexp]->eval(m[i]->u, m[i]->v);
	chi2 += pow(val - m[i]->mag_cat, 2.0) / (pow(m[i]->err, 2.0) + pow(m[i]->err_cat, 2.0));
	mag2 += pow(val - m[i]->mag_cat, 2.0);
	num++;
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;
	double val = s[i]->mag - 2.5 * log10(fexp[s[i]->iexp] * fchip[s[i]->ichip]);
	val += ffpSet[s[i]->iexp]->eval(s[i]->u, s[i]->v);
	chi2 += pow((val - s[i]->mag0)/s[i]->err, 2.0);
	mag2 += pow((val - s[i]->mag0), 2.0);
	num++;
    }

    if (mag)
	return mag2 / num;
    else
	return chi2 / num;
}

void flagObj_rel(std::vector<Obs::Ptr> &m,
		 std::vector<Obs::Ptr> &s,
		 double e2,
		 std::map<int, float>& fexp,
		 std::map<int, float>& fchip,
		 //FluxFitParams::Ptr p)
		 FfpSet &ffpSet)
{
    int nMobs = m.size();
    int nSobs = s.size();

    int nreject = 0;
    for (int i = 0; i < nMobs; i++) {
	if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 || m[i]->err == -9999) continue;
	double val = m[i]->mag - 2.5 * log10(fexp[m[i]->iexp] * fchip[m[i]->ichip]);
	val += ffpSet[m[i]->iexp]->eval(m[i]->u, m[i]->v);
	//double r2 = pow((val - m[i]->mag0)/m[i]->err, 2.0);
	double r2 = pow((val - m[i]->mag0), 2.0);
	if (r2 > e2) {
	    m[i]->good = false;
	    nreject++;
	}
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;
	double val = s[i]->mag - 2.5 * log10(fexp[s[i]->iexp] * fchip[s[i]->ichip]);
	val += ffpSet[s[i]->iexp]->eval(s[i]->u, s[i]->v);
	//double r2 = pow((val - s[i]->mag0)/s[i]->err, 2.0);
	double r2 = pow((val - s[i]->mag0), 2.0);
	if (r2 > e2) {
	    s[i]->good = false;
	    nreject++;
	}
    }

    printf("nreject: %d\n", nreject);
}

void flagObj_abs(std::vector<Obs::Ptr> &m,
		 std::vector<Obs::Ptr> &s,
		 double e2,
		 std::map<int, float>& fexp,
		 std::map<int, float>& fchip,
		 //FluxFitParams::Ptr p)
		 FfpSet &ffpSet)
{
    int nMobs = m.size();
    int nSobs = s.size();

    int nreject = 0;
    for (int i = 0; i < nMobs; i++) {
	if (m[i]->jstar == -1 || !m[i]->good || m[i]->mag == -9999 ||
	    m[i]->err == -9999 || m[i]->mag_cat == -9999) continue;
	double val = m[i]->mag - 2.5 * log10(fexp[m[i]->iexp] * fchip[m[i]->ichip]);
	//val += p->eval(m[i]->u, m[i]->v);
	val += ffpSet[m[i]->iexp]->eval(m[i]->u, m[i]->v);
	//double r2 = pow(val - m[i]->mag_cat, 2.0) / (pow(m[i]->err, 2.0) + pow(m[i]->err_cat, 2.0));
	//if (r2 > e2) {
	double r2 = pow(val - m[i]->mag_cat, 2.0);
	if (r2 > e2 + 9.0 * pow(m[i]->err_cat, 2.0)) {
	    m[i]->good = false;
	    nreject++;
	}
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999 || s[i]->err == -9999) continue;
	double val = s[i]->mag - 2.5 * log10(fexp[s[i]->iexp] * fchip[s[i]->ichip]);
	//val += p->eval(s[i]->u, s[i]->v);
	val += ffpSet[s[i]->iexp]->eval(s[i]->u, s[i]->v);
	//double r2 = pow((val - s[i]->mag0)/s[i]->err, 2.0);
	double r2 = pow((val - s[i]->mag0), 2.0);
	if (r2 > e2) {
	    s[i]->good = false;
	    nreject++;
	}
    }

    printf("nreject: %d\n", nreject);
}

void fluxFitRelative(ObsVec& matchVec,
		     int nmatch,
		     ObsVec& sourceVec,
		     int nsource,
		     WcsDic& wcsDic,
		     CcdSet& ccdSet,
		     std::map<int, float>& fexp,
		     std::map<int, float>& fchip,
		     FfpSet &ffpSet,
		     bool solveCcd,
		     bool common) {

    int nexp = wcsDic.size();
    int nchip = ccdSet.size();

    for (int k = 0; k < 3; k++) {
	Eigen::VectorXd fsol;
	for (WcsDic::iterator it = wcsDic.begin(); it != wcsDic.end(); it++) {
	    ffpSet[it->first]->coeff[0] = 0.0;
	}
	if (common) {
	    fsol = fluxFit_rel1(matchVec, nmatch, sourceVec, nsource, nexp, nchip, ffpSet, solveCcd);
	} else {
	    fsol = fluxFit_rel(matchVec, nmatch, sourceVec, nsource, nexp, nchip, ffpSet, solveCcd);
	}
	int i = 0;
	for (WcsDic::iterator it = wcsDic.begin();
	     it != wcsDic.end(); it++, i++) {
	    //fexp.insert(std::map<int, float>::value_type(it->first, pow(10., -0.4*fsol(i))));
	    fexp[it->first] = pow(10., -0.4*fsol(i));
	    // Set the flux correction value at the origin of focal plane coordinate to be zero
	    double f0 = ffpSet[it->first]->eval(0,0);
	    ffpSet[it->first]->coeff[0] = -f0;
	    fexp[it->first] = pow(10., -0.4*(fsol(i)+f0));
	}
	if (solveCcd) {
	    for (CcdSet::iterator it = ccdSet.begin();
		 it != ccdSet.end(); it++, i++) {
		//fchip.insert(std::map<int, float>::value_type(it->first, pow(10., -0.4*fsol(i))));
		fchip[it->first] = pow(10., -0.4*fsol(i));
	    }
	} else {
	    for (CcdSet::iterator it = ccdSet.begin();
		 it != ccdSet.end(); it++, i++) {
		//fchip.insert(std::map<int, float>::value_type(it->first, 1.0));
		fchip[it->first] = 1.0;
	    }
	}
	double chi2f = calcChi2_rel(matchVec, sourceVec, fexp, fchip, ffpSet);
	printf("chi2f: %e\n", chi2f);
	double e2f = calcChi2_rel(matchVec, sourceVec, fexp, fchip, ffpSet, true);
	printf("err: %f (mag)\n", sqrt(e2f));
	//flagObj_rel(matchVec, sourceVec, 9.0, fexp, fchip, ffp);
	if (k < 2)
	    flagObj_rel(matchVec, sourceVec, 9.0*e2f, fexp, fchip, ffpSet);
    }

    printf("FFP:   ");
    for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++) {
	printf(" %8d", it->first);
    }
    printf("\n");
    printf("FFP:   ");
    for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++) {
	printf(" %8.5f", ffpSet[it->first]->eval(0,0));
    }
    printf("\n");
    FluxFitParams::Ptr ffp = ffpSet[ffpSet.begin()->first];
    for (int i = 0; i < ffp->ncoeff; i++) {
	printf("FFP: %2d", i);
	for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++) {
	    printf(" %8.5f", ffpSet[it->first]->coeff[i]);
	}
	printf("\n");
    }

    for (typename std::map<int, float>::iterator it = fchip.begin(); it != fchip.end(); ++it) {
        std::cout << "CCD " << it->first << ": " << it->second << std::endl;
    }
}

void fluxFitAbsolute(ObsVec& matchVec,
		     int nmatch,
		     ObsVec& sourceVec,
		     int nsource,
		     WcsDic& wcsDic,
		     CcdSet& ccdSet,
		     std::map<int, float>& fexp,
		     std::map<int, float>& fchip,
		     FfpSet &ffpSet,
		     bool solveCcd,
                     bool common) {

    int nexp = wcsDic.size();
    int nchip = ccdSet.size();

    for (int k = 0; k < 3; k++) {
	Eigen::VectorXd fsol;
	for (WcsDic::iterator it = wcsDic.begin(); it != wcsDic.end(); it++) {
	    ffpSet[it->first]->coeff[0] = 0.0;
	}
	if (common) {
	    fsol = fluxFit_abs1(matchVec, nmatch, sourceVec, nsource, nexp, nchip, ffpSet, solveCcd);
	} else {
	    fsol = fluxFit_abs(matchVec, nmatch, sourceVec, nsource, nexp, nchip, ffpSet, solveCcd);
	}
	int i = 0;
	for (WcsDic::iterator it = wcsDic.begin();
	     it != wcsDic.end(); it++, i++) {
	    //fexp.insert(std::map<int, float>::value_type(it->first, pow(10., -0.4*fsol(i))));
	    fexp[it->first] = pow(10., -0.4*fsol(i));
	    // Set the flux correction value at the origin of focal plane coordinate to be zero
	    double f0 = ffpSet[it->first]->eval(0,0);
	    ffpSet[it->first]->coeff[0] = -f0;
	    fexp[it->first] = pow(10., -0.4*(fsol(i)+f0));
	}
	if (solveCcd) {
	    for (CcdSet::iterator it = ccdSet.begin();
		 it != ccdSet.end(); it++, i++) {
		//fchip.insert(std::map<int, float>::value_type(it->first, pow(10., -0.4*fsol(i))));
		fchip[it->first] = pow(10., -0.4*fsol(i));
	    }
	} else {
	    for (CcdSet::iterator it = ccdSet.begin();
		 it != ccdSet.end(); it++, i++) {
		//fchip.insert(std::map<int, float>::value_type(it->first, 1.0));
		fchip[it->first] = 1.0;
	    }
	}
	double chi2f = calcChi2_abs(matchVec, sourceVec, fexp, fchip, ffpSet);
	printf("chi2f: %e\n", chi2f);
	double e2f = calcChi2_abs(matchVec, sourceVec, fexp, fchip, ffpSet, true);
	printf("err: %f (mag)\n", sqrt(e2f));
	//flagObj_abs(matchVec, sourceVec, 9.0, fexp, fchip, ffp);
	if (k < 2)
	    flagObj_abs(matchVec, sourceVec, 9.0*e2f, fexp, fchip, ffpSet);
    }

    printf("FFP:   ");
    for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++) {
	printf(" %8d", it->first);
    }
    printf("\n");
    printf("FFP:   ");
    for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++) {
	printf(" %8.5f", ffpSet[it->first]->eval(0,0));
    }
    printf("\n");
    FluxFitParams::Ptr ffp = ffpSet[ffpSet.begin()->first];
    for (int i = 0; i < ffp->ncoeff; i++) {
	printf("FFP: %2d", i);
	for (FfpSet::iterator it = ffpSet.begin(); it != ffpSet.end(); it++) {
	    printf(" %8.5f", ffpSet[it->first]->coeff[i]);
	}
	printf("\n");
    }

    for (typename std::map<int, float>::iterator it = fchip.begin(); it != fchip.end(); ++it) {
        std::cout << "CCD " << it->first << ": " << it->second << std::endl;
    }
}

void lsst::meas::mosaic::fluxFit(bool absolute,
				 bool common,
				 ObsVec& matchVec,
				 int nmatch,
				 ObsVec& sourceVec,
				 int nsource,
				 WcsDic& wcsDic,
				 CcdSet& ccdSet,
				 std::map<int, float>& fexp,
				 std::map<int, float>& fchip,
				 FfpSet &ffpSet,
				 bool solveCcd) {
    printf("fluxFit ...\n");
    if (absolute) {
	fluxFitAbsolute(matchVec, nmatch, sourceVec, nsource, wcsDic, ccdSet, fexp, fchip, ffpSet, solveCcd, common);
    } else {
	fluxFitRelative(matchVec, nmatch, sourceVec, nsource, wcsDic, ccdSet, fexp, fchip, ffpSet, solveCcd, common);
    }
}

FluxFitParams::Ptr
lsst::meas::mosaic::convertFluxFitParams(FluxFitParams::Ptr& ffp, PTR(lsst::afw::cameraGeom::Detector)& ccd, double x0, double y0)
{
    FluxFitParams::Ptr newP = FluxFitParams::Ptr(new FluxFitParams(ffp->order, ffp->chebyshev));
    newP->u_max = 1.0;
    newP->v_max = 1.0;

    int *xorder = ffp->xorder;
    int *yorder = ffp->yorder;

    lsst::afw::cameraGeom::Orientation ori = ccd->getOrientation();
    double cosYaw = ori.getCosYaw();
    double sinYaw = ori.getSinYaw();

    // u = cc * u' - ss * v'
    // v = ss * u' + cc * v'
    // u^i * v^j = (cc * u' - ss * v')^i * (ss * u' + cc * v')^j
    //           = \Sigma (i, n) * (cc * u')^n * (-ss * v')^(i-n) *
    //             \Sigma (j, m) * (ss * u')^m * ( cc * v')^(j-m)
    for (int k = 0; k < ffp->ncoeff; k++) {
	for (int n = 0; n <= xorder[k]; n++) {
	    for (int m = 0; m <= yorder[k]; m++) {
		int i = n + m;
		int j = xorder[k] + yorder[k] - n - m;
		int l = newP->getIndex(i, j);
		double C =  binomial(xorder[k], n) *
		            binomial(yorder[k], m) *
		            pow(cosYaw, n) * pow(-sinYaw, xorder[k]-n) *
		            pow(sinYaw, m) * pow( cosYaw, yorder[k]-m)
		            / pow(ffp->u_max, xorder[k])
		            / pow(ffp->v_max, yorder[k]);
		newP->coeff[l] += ffp->coeff[k] * C;
	    }
	}
    }

    afw::geom::Extent2D off = ccd->getCenter().getPixels(ccd->getPixelSize()) - ccd->getCenterPixel();
    newP->x0 =  (off[0] + x0) * cosYaw + (off[1] + y0) * sinYaw;
    newP->y0 = -(off[0] + x0) * sinYaw + (off[1] + y0) * cosYaw;

    return newP;
}

lsst::daf::base::PropertySet::Ptr
lsst::meas::mosaic::metadataFromFluxFitParams(FluxFitParams::Ptr& ffp)
{
    lsst::daf::base::PropertySet::Ptr metadata =
       lsst::daf::base::PropertySet::Ptr(new lsst::daf::base::PropertySet());

    metadata->set("ORDER", ffp->order);
    metadata->set("ABSOLUTE", ffp->absolute);
    metadata->set("CHEBYSHEV", ffp->chebyshev);
    metadata->set("NCOEFF", ffp->ncoeff);
    metadata->set("U_MAX", ffp->u_max);
    metadata->set("V_MAX", ffp->v_max);
    metadata->set("X0", ffp->x0);
    metadata->set("Y0", ffp->y0);

    for (int k = 0; k < ffp->ncoeff; k++) {
       std::string label = boost::str(boost::format("C_%d_%d") % ffp->xorder[k] % ffp->yorder[k]);
       metadata->set(label, ffp->coeff[k]);
    }

    return metadata;
}

lsst::afw::image::Image<float>::Ptr
lsst::meas::mosaic::getFCorImg(FluxFitParams::Ptr& p,
			      PTR(lsst::afw::cameraGeom::Detector)& ccd,
			      Coeff::Ptr& coeff)
{
    int width  = ccd->getAllPixels(true).getWidth();
    int height = ccd->getAllPixels(true).getHeight();

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

        afw::geom::Point2D uv 
            = ccd->getPositionFromPixel(afw::geom::Point2D(x, y)).getPixels(ccd->getPixelSize())
            + afw::geom::Extent2D(coeff->x0, coeff->y0);
	    double val0 = p->eval(uv.getX(), uv.getY());
        uv = ccd->getPositionFromPixel(afw::geom::Point2D(xend, y)).getPixels(ccd->getPixelSize())
            + afw::geom::Extent2D(coeff->x0, coeff->y0);
	    double val1 = p->eval(uv.getX(), uv.getY());

	    for (int i = 0; i < interval; i++) {
		vals(x+i) = val0 + (val1 - val0) / interval * i;
	    }
	}

	lsst::afw::image::Image<float>::x_iterator begin = img->row_begin(y);
	lsst::afw::image::Image<float>::x_iterator end   = img->row_end(y);

	for (lsst::afw::image::Image<float>::x_iterator ptr = begin; ptr != end; ptr++) {

	    int x = ptr - begin;

	    *ptr = pow(10., -0.4*vals(x));
	}
    }

    return img;
}

lsst::afw::image::Image<float>::Ptr
lsst::meas::mosaic::getFCorImg(FluxFitParams::Ptr& p, int width, int height)
{
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
	    double val0 = p->eval(u, v);
	    u = xend;
	    v = y;
	    double val1 = p->eval(u, v);
	    for (int i = 0; i < interval; i++) {
		vals(x+i) = val0 + (val1 - val0) / interval * i;
	    }
	}

	lsst::afw::image::Image<float>::x_iterator begin = img->row_begin(y);
	lsst::afw::image::Image<float>::x_iterator end   = img->row_end(y);

	for (lsst::afw::image::Image<float>::x_iterator ptr = begin; ptr != end; ptr++) {

	    int x = ptr - begin;

	    *ptr = pow(10., -0.4*vals(x));
	}
    }

    return img;
}

lsst::afw::image::Image<float>::Ptr
lsst::meas::mosaic::getFCorImg(FluxFitParams::Ptr& p,
			      PTR(lsst::afw::cameraGeom::Detector)& ccd)
{
    int width  = ccd->getAllPixels(true).getWidth();
    int height = ccd->getAllPixels(true).getHeight();

    return getFCorImg(p, width, height);
}
