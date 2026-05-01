/**
 * @file legendre_rule.h
 * @brief Public API for Gauss-Legendre quadrature rule generation.
 *
 * Provides a single function that computes the nodes and weights of
 * an n-point Gauss-Legendre quadrature rule on an arbitrary interval
 * [a, b].
 */
 
#ifndef LEGENDRE_RULE_H
#define LEGENDRE_RULE_H
 
#ifdef __cplusplus
extern "C" {
#endif
 
/**
 * @brief Compute an n-point Gauss-Legendre quadrature rule on [a, b].
 *
 * After the call, the rule approximates
 * @f[
 *     \int_a^b f(t)\,dt \;\approx\; \sum_{i=0}^{n-1} w_i \, f(x_i).
 * @f]
 *
 * The caller is responsible for allocating @p x and @p w with at
 * least @p n elements each.
 *
 * @param[in]  n  Number of quadrature points (must be >= 1).
 * @param[in]  a  Left endpoint of the interval.
 * @param[in]  b  Right endpoint of the interval.
 * @param[out] x  Array of length @p n; receives the abscissas (nodes).
 * @param[out] w  Array of length @p n; receives the weights.
 *
 * @note Nodes are returned in ascending order. Weights are positive
 *       and sum to (b - a).
 */
void legendre_rule(int n, double a, double b, double x[], double w[]);
 
#ifdef __cplusplus
}
#endif
 
#endif /* LEGENDRE_RULE_H */