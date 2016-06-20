/**
 * @file Predictor.h
 * @brief Dynamics Predictor
 * @author Paul Drews
 */

#pragma once

#include <gtsam_unstable/base/dllexport.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/nonlinear/Expression.h>
#include <math.h>

#define X_DIM 5
#define U_DIM 2
#define THETA_DIM 130

namespace gtsam {
//Convenience typedef
typedef Expression<double> double_;

/**
 * Create an expression out of a linear function f:T->A with (constant) Jacobian dTdA
 * Specialize for stacking scalars to vector
 * T is the return vector type, A is the scalar type of the expression
 * Only works for scalar A at the moment
 */
template <typename T, typename A>
Expression<T> StackExpression2(
    const Expression<A>& expression1, const Expression<A>& expression2);

/**
 * Create an expression out of a linear function f:T->A with (constant) Jacobian dTdA
 * Specialize for stacking scalars to vector
 * T is the return vector type, A is the scalar type of the expression
 */
template <typename T, typename A>
Expression<T> StackExpression3(
    const Expression<A>& expression1, const Expression<A>& expression2, const Expression<A>& expression3);

/**
 * Create a cosine factor
 */
double_ cosExpression(const double_& expression);

/**
 * Create a cosine factor
 */
double_ sinExpression(const double_& expression);

double_ arctanExpression(const double_& expression);

double_ tanh_(const double_& expression);

double_ sqrt_(const double_& expression);

double_ Mul_(
    const double_& expression1, const double_& expression2);

double_ Div_(
    const double_& expression1, const double_& expression2);

double_ gammaExpression(const double_& mu,
        const double_& Fz, const double_& U2);

double_ FyrExpression(const double_& mu,
    const double_& Fz, const double_& gamma,
    const double_& cAlpha, const double_& alpha);

double_ FyfExpression(const double_& mu,
    const double_& Fz,
    const double_& cAlpha, const double_& alpha);

double_ betadotExpression(const double_& mu,
    const double_& m, const double_& Vx,
    const double_& thetad, const double_& Fyf,
    const double_& Fyr);

double_ VxdExpression(const double_& U2,
    const double_& Fyf, const double_& U1,
    const double_& m, const double_& Vx,
    const double_& beta, const double_& thetad);

double_ thetaddExpression(const double_& a,
    const double_& Fyf, const double_& b,
    const double_& Fyr, const double_& Iz);

double_ alphafExpression(const double_& beta,
    const double_& a, const double_& Vx,
    const double_& thetad, const double_& delta);

double_ alpharExpression(const double_& beta,
    const double_& a, const double_& Vx,
    const double_& thetad);

double_ alphaslExpression(const double_& mu,
    const double_& Fz, const double_& gamma,
    const double_& cAlpha);


} // \namespace gtsam
