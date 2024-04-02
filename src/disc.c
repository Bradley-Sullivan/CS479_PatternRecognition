#include "disc.h"

/**
 * @brief Simple Euclidean distance classifier. Optimum
 * for classes with equal, diagonal covariance matrices
 * which all share the same variance and priors. Defines
 * a linear decision boundary.
 * 
 * @param x - Feature vector
 * @param g - Class/Category distribution
 * @return double - (-1) * dist_from_boundary
 */
double euclid_disc(VEC *x, gauss_t *g) {
    VEC *diff = v_sub(x, g->mu, VNULL);
    return -1 * in_prod(diff, diff);
}

/**
 * @brief Optimum classifier for m.v. Gaussian classes
 * with equal, diagonal covariance matrices w/ same variance.
 * Simplifies to Euclidean distance classifier for classes
 * with equal priors. Linear decision boundary.
 * 
 * @param x - Feature vector
 * @param g - Class/Category distribution
 * @return double - (-1) * dist_from_boundary
 */
double case1_disc(VEC *x, gauss_t *g) {
    
    // g(x) = wT.x + w0
    
    VEC *w;
    double w0, ss, n;

    ss = g->sigma->me[0][0] * g->sigma->me[0][0];

    w = sv_mlt(1 / ss, g->mu, VNULL);
    w0 = (-1 / (2 * ss)) * in_prod(g->mu, g->mu) + log(g->prior);

    n = in_prod(w, x) + w0;

    return n;
}

/**
 * @brief Optimum classifier for m.v. Gaussian classes,
 * each with arbitrary covariance matrices. Hyperquadric
 * decision boundary.
 * 
 * @param x - Feature vector
 * @param g - Class/Category distribution
 * @return double - (-1) * dist_from_boundary
 */
double case3_disc(VEC *x, gauss_t *g) {

    // g(x) = xT.W.x + wT.x + wi0

    double n, w0, det;
    MAT *W, *inv;
    VEC *w, *a;

    inv = m_inverse(g->sigma, MNULL);
    W = sm_mlt(-0.5, inv, MNULL);

    w = mv_mlt(inv, g->mu, VNULL);

    w0 = -0.5 * in_prod(g->mu, w);
    w0 +=  -0.5 * log(m_det(g->sigma));
    w0 += log(g->prior);

    a = mv_mlt(W, x, VNULL);

    n = in_prod(x, a) + in_prod(w, x) + w0;

    m_free(inv); m_free(W);
    v_free(w); v_free(a);

    return n;
}