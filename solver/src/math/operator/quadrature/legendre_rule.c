/**
 * @file legendre_rule.c
 * @brief Pure-function implementation of the Glaser-Liu-Rokhlin
 *        Gauss-Legendre quadrature rule generator.
 *
 * Stripped of file I/O, timing, and the original main(). The numerical
 * core is unchanged from the reference implementation and exposes a
 * single public entry point: legendre_rule().
 *
 * @par Reference
 *      A. Glaser, X. Liu, V. Rokhlin,
 *      "A fast algorithm for the calculation of the roots of special
 *      functions," SIAM Journal on Scientific Computing, 29(4),
 *      pp. 1420-1438, 2007.
 *
 * @par Original authors
 *      Nick Hale and John Burkardt. Distributed under the GNU LGPL.
 * 
 * @link https://people.math.sc.edu/Burkardt/c_src/legendre_rule_fast/legendre_rule_fast.html
 *       https://people.math.sc.edu/Burkardt/c_src/legendre_rule_fast/legendre_rule_fast.c
 */

#include <stdlib.h>
#include <math.h>

#include "math/operator/quadrature/legendre_rule.h"

/* ---------- internal helpers ---------- */

static void   legendre_compute_glr  (int n, double x[], double w[]);
static void   legendre_compute_glr0 (int n, double *p, double *pp);
static void   legendre_compute_glr1 (int n, double *roots, double *ders);
static void   legendre_compute_glr2 (double p, int n, double *roots, double *ders);
static void   rescale               (double a, double b, int n, double x[], double w[]);
static double rk2_leg               (double t1, double t2, double x, int n);
static double ts_mult               (double *u, double h, int n);

/* ---------- public entry point ---------- */

/**
 * @brief Compute an n-point Gauss-Legendre quadrature rule on [a, b].
 *
 * Internally this calls legendre_compute_glr() to obtain the rule on
 * the canonical interval [-1, 1], then rescales nodes and weights to
 * [a, b] via rescale().
 *
 * @param[in]  n  Number of quadrature points (must be >= 1).
 * @param[in]  a  Left endpoint of the interval.
 * @param[in]  b  Right endpoint of the interval.
 * @param[out] x  Array of length @p n; receives the abscissas (nodes).
 * @param[out] w  Array of length @p n; receives the weights.
 *
 * @see legendre_compute_glr()
 * @see rescale()
 */
void legendre_rule(int n, double a, double b, double x[], double w[])
{
    if (n < 1 || x == NULL || w == NULL) {
        return;
    }
    legendre_compute_glr(n, x, w);
    rescale(a, b, n, x, w);
}

/* ---------- algorithm (unchanged numerics) ---------- */

/**
 * @brief Compute Gauss-Legendre nodes and weights on [-1, 1] via GLR.
 *
 * Implements the Glaser-Liu-Rokhlin fast algorithm. For odd @p n the
 * polynomial value at zero is itself a root; for even @p n the first
 * positive root is found via Newton iteration on a Taylor expansion.
 * The remaining roots are then propagated by sequential Taylor
 * expansion along the chain of roots, and the weights are recovered
 * from the derivatives at the nodes.
 *
 * @param[in]  n  Order of the rule.
 * @param[out] x  Array of length @p n; receives the abscissas.
 * @param[out] w  Array of length @p n; receives the weights.
 */
static void legendre_compute_glr(int n, double x[], double w[])
{
    int i;
    double p, pp, w_sum;

    legendre_compute_glr0(n, &p, &pp);

    /* Either zero is a root, or we must locate the first positive root. */
    if (n % 2 == 1) {
        x[(n - 1) / 2] = p;
        w[(n - 1) / 2] = pp;
    } else {
        legendre_compute_glr2(p, n, &x[n / 2], &w[n / 2]);
    }

    legendre_compute_glr1(n, x, w);

    /* Convert derivatives stored in w[] into actual quadrature weights. */
    for (i = 0; i < n; i++)
    {
        w[i] = 2.0 / (1.0 - x[i]) / (1.0 + x[i]) / w[i] / w[i];
    }
    w_sum = 0.0;
    for (i = 0; i < n; i++)
    {
        w_sum += w[i];
    }
    for (i = 0; i < n; i++)
    {
        w[i] = 2.0 * w[i] / w_sum;
    }
}

/**
 * @brief Compute P_n(0) and P_n'(0) by upward recursion.
 *
 * Uses the standard three-term recurrence for the Legendre polynomial
 * and its derivative, evaluated at t = 0.
 *
 * @param[in]  n   Order of the Legendre polynomial.
 * @param[out] p   Receives P_n(0).
 * @param[out] pp  Receives P_n'(0).
 */
static void legendre_compute_glr0(int n, double *p, double *pp)
{
    double dk;
    int k;
    double pm1, pm2, ppm1, ppm2;

    pm2 = 0.0;
    pm1 = 1.0;
    ppm2 = 0.0;
    ppm1 = 0.0;

    for (k = 0; k < n; k++)
    {
        dk = (double) k;
        *p  = -dk * pm2 / (dk + 1.0);
        *pp = ((2.0 * dk + 1.0) * pm1 - dk * ppm2) / (dk + 1.0);
        pm2 = pm1;
        pm1 = *p;
        ppm2 = ppm1;
        ppm1 = *pp;
    }
}

/**
 * @brief Propagate roots and derivatives via Taylor expansion.
 *
 * Given a starting estimate at index @c n2 (where @c n2 = (n-1)/2 for
 * odd @p n, or @c n/2 for even @p n), this routine builds the
 * remaining roots in the upper half by Taylor-series stepping plus
 * Newton refinement, then mirrors them into the lower half by the
 * symmetry of P_n.
 *
 * @param[in]     n     Order of the Legendre polynomial.
 * @param[in,out] x     On input, contains a starting root in entry
 *                      @c n2. On output, contains all @p n roots.
 * @param[in,out] ders  On input, contains the derivative at the
 *                      starting root. On output, contains the
 *                      derivatives at all @p n roots.
 */
static void legendre_compute_glr1(int n, double *x, double *ders)
{
    double dk, dn, h;
    int j, k, l;
    int m = 30;       /**< number of terms in the Taylor expansion */
    int n2, s;
    const double pi = 3.141592653589793;
    double *u, *up;
    double xp;

    if (n % 2 == 1) {
        n2 = (n - 1) / 2;
        s = 1;
    } else {
        n2 = n / 2;
        s = 0;
    }

    u  = (double *) malloc((m + 2) * sizeof(double));
    up = (double *) malloc((m + 1) * sizeof(double));

    dn = (double) n;

    for (j = n2; j < n - 1; j++)
    {
        xp = x[j];

        h = rk2_leg(pi / 2.0, -pi / 2.0, xp, n) - xp;

        u[0] = 0.0;
        u[1] = 0.0;
        u[2] = ders[j];

        up[0] = 0.0;
        up[1] = u[2];

        /* Build the Taylor coefficients of P_n around xp. */
        for (k = 0; k <= m - 2; k++)
        {
            dk = (double) k;
            u[k + 3] = (
                2.0 * xp * (dk + 1.0) * u[k + 2]
                + (dk * (dk + 1.0) - dn * (dn + 1.0)) * u[k + 1] / (dk + 1.0)
            ) / (1.0 - xp) / (1.0 + xp) / (dk + 2.0);

            up[k + 2] = (dk + 2.0) * u[k + 3];
        }

        /* Newton iteration on the Taylor series to refine the next root. */
        for (l = 0; l < 5; l++)
        {
            h = h - ts_mult(u, h, m) / ts_mult(up, h, m - 1);
        }

        x[j + 1]    = xp + h;
        ders[j + 1] = ts_mult(up, h, m - 1);
    }

    free(u);
    free(up);

    /* Mirror the upper-half roots to the lower half. */
    for (k = 0; k < n2 + s; k++)
    {
        x[k]    = -x[n - k - 1];
        ders[k] =  ders[n - k - 1];
    }
}

/**
 * @brief Find the first positive root of P_n when @p n is even.
 *
 * Combines a Runge-Kutta integration of an ODE for the root location
 * (rk2_leg()) with Newton refinement on a Taylor expansion of P_n.
 *
 * @param[in]  pn0  Value of P_n at zero.
 * @param[in]  n    Order of the Legendre polynomial.
 * @param[out] x1   Receives the first positive root.
 * @param[out] d1   Receives the derivative P_n'(x1).
 *
 * @note Only called when @p n is even; for odd @p n, zero is itself
 *       a root and this routine is bypassed.
 */
static void legendre_compute_glr2(double pn0, int n, double *x1, double *d1)
{
    double dk, dn;
    int k, l;
    int m = 30;       /**< number of terms in the Taylor expansion */
    const double pi = 3.141592653589793;
    double t;
    double *u, *up;

    t = 0.0;
    *x1 = rk2_leg(t, -pi / 2.0, 0.0, n);

    u  = (double *) malloc((m + 2) * sizeof(double));
    up = (double *) malloc((m + 1) * sizeof(double));

    dn = (double) n;

    u[0] = 0.0;
    u[1] = pn0;
    up[0] = 0.0;

    for (k = 0; k <= m - 2; k += 2)
    {
        dk = (double) k;

        u[k + 2] = 0.0;
        u[k + 3] = (dk * (dk + 1.0) - dn * (dn + 1.0)) * u[k + 1]
                   / (dk + 1.0) / (dk + 2.0);

        up[k + 1] = 0.0;
        up[k + 2] = (dk + 2.0) * u[k + 3];
    }

    /* Newton refinement of the root estimate. */
    for (l = 0; l < 5; l++)
    {
        *x1 = *x1 - ts_mult(u, *x1, m) / ts_mult(up, *x1, m - 1);
    }
    *d1 = ts_mult(up, *x1, m - 1);

    free(u);
    free(up);
}

/**
 * @brief Affinely rescale a quadrature rule from [-1, 1] to [a, b].
 *
 * Applies the linear map x' = ((a + b) + (b - a) x) / 2 to the nodes
 * and w' = (b - a) w / 2 to the weights, in place.
 *
 * @param[in]     a  New left endpoint.
 * @param[in]     b  New right endpoint.
 * @param[in]     n  Number of nodes.
 * @param[in,out] x  Nodes; rescaled in place.
 * @param[in,out] w  Weights; rescaled in place.
 */
static void rescale(double a, double b, int n, double x[], double w[])
{
    int i;

    for (i = 0; i < n; i++) 
    {
        x[i] = ((a + b) + (b - a) * x[i]) / 2.0;
    }
    for (i = 0; i < n; i++) 
    {
        w[i] = (b - a) * w[i] / 2.0;
    }
}

/**
 * @brief Advance an ODE for a Legendre root via two-stage Runge-Kutta.
 *
 * Integrates the ODE that tracks the root location of P_n along an
 * angular parameterization. Used to obtain a high-quality starting
 * estimate before Newton refinement.
 *
 * @param[in] t1  Start of the integration interval.
 * @param[in] t2  End of the integration interval.
 * @param[in] x   Value of x at @p t1.
 * @param[in] n   Order of the Legendre polynomial.
 * @return        Value of x at @p t2.
 */
static double rk2_leg(double t1, double t2, double x, int n)
{
    double f, h, k1, k2, snn1, t;
    int j;
    int m = 10;       /**< number of RK steps */

    h = (t2 - t1) / (double) m;
    snn1 = sqrt((double)(n * (n + 1)));

    t = t1;

    for (j = 0; j < m; j++)
    {
        f = (1.0 - x) * (1.0 + x);
        k1 = -h * f / (snn1 * sqrt(f) - 0.5 * x * sin(2.0 * t));
        x = x + k1;

        t = t + h;

        f = (1.0 - x) * (1.0 + x);
        k2 = -h * f / (snn1 * sqrt(f) - 0.5 * x * sin(2.0 * t));
        x = x + 0.5 * (k2 - k1);
    }
    return x;
}

/**
 * @brief Evaluate a polynomial in monomial form, ignoring @c u[0].
 *
 * Computes
 * @f[
 *     u_1 + u_2\, h + u_3\, h^2 + \dots + u_n\, h^{n-1}.
 * @f]
 *
 * @param[in] u  Coefficient array (@c u[0] is unused).
 * @param[in] h  Argument.
 * @param[in] n  Number of terms to sum.
 * @return       The polynomial value.
 */
static double ts_mult(double *u, double h, int n)
{
    double hk, ts;
    int k;

    ts = 0.0;
    hk = 1.0;
    for (k = 1; k <= n; k++)
    {
        ts += u[k] * hk;
        hk *= h;
    }
    return ts;
}