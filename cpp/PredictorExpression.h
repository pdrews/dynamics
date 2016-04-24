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
    const Expression<A>& expression1, const Expression<A>& expression2) {
  typename Expression<T>::template BinaryFunction<A,A>::type g =
      [=](const A& value1, const A& value2, typename MakeOptionalJacobian<T, A>::type H1, typename MakeOptionalJacobian<T, A>::type H2) {
        if (H1)
          *H1 << 1,0;
        if (H2)
          *H2 << 0,1;
        T ret;
        ret << value1, value2;
        std::cout << ret << std::endl << value1 << std::endl << value2 << std::endl;
        return ret;
      };
  return Expression<T>(g, expression1, expression2);
}

/**
 * Create an expression out of a linear function f:T->A with (constant) Jacobian dTdA
 * Specialize for stacking scalars to vector
 * T is the return vector type, A is the scalar type of the expression
 */
template <typename T, typename A>
Expression<T> StackExpression3(
    const Expression<A>& expression1, const Expression<A>& expression2, const Expression<A>& expression3) {
  typename Expression<T>::template TernaryFunction<A,A,A>::type g =
      [=](const A& value1, const A& value2, const A& value3,
          typename MakeOptionalJacobian<T, A>::type H1, typename MakeOptionalJacobian<T, A>::type H2, typename MakeOptionalJacobian<T, A>::type H3) {
        if (H1)
          *H1 << 1,0,0;
        if (H2)
          *H2 << 0,1,0;
        if (H3)
          *H3 << 0,0,1;
        T ret;
        ret << value1, value2, value3;
        return ret;
      };
  return Expression<T>(g, expression1, expression2, expression3);
}

/**
 * Create a cosine factor
 */
double_ cosExpression(const double_& expression) {
  // Use lambda to calculate the Jacobian
  typename double_::template UnaryFunction<double>::type g =
      [=](const double& value, OptionalJacobian<1,1> H) {
        if (H)
          *H << -sin(value);
        return cos(value);
      };
  return double_(g, expression);
}

/**
 * Create a cosine factor
 */
double_ sinExpression(const double_& expression) {
  // Use lambda to calculate the Jacobian
  typename double_::template UnaryFunction<double>::type g =
      [=](const double& value, OptionalJacobian<1,1> H) {
        if (H)
          *H << cos(value);
        return sin(value);
      };
  return double_(g, expression);
}

double_ arctanExpression(const double_& expression) {
  // Use lambda to calculate the Jacobian
  typename double_::template UnaryFunction<double>::type g =
      [=](const double& value, OptionalJacobian<1,1> H) {
        if (H)
          *H << 1/((value*value)+1);
        return atan(value);
      };
  return double_(g, expression);
}

double_ tanhExpression(const double_& expression) {
  // Use lambda to calculate the Jacobian
  typename double_::template UnaryFunction<double>::type g =
      [=](const double& value, OptionalJacobian<1,1> H) {
        if (H)
          *H << pow(1/(cosh(value)),2);
        return tanh(value);
      };
  return double_(g, expression);
}

double_ sqrtExpression(const double_& expression) {
  // Use lambda to calculate the Jacobian
  typename double_::template UnaryFunction<double>::type g =
      [=](const double& value, OptionalJacobian<1,1> H) {
        if (H) {
          assert(value > 0);
          *H << 1/(2*sqrt(value));
        }
        return sqrt(value);
      };
  return double_(g, expression);
}

double_ MulExpression(
    const double_& expression1, const double_& expression2) {
  typename double_::template BinaryFunction<double,double>::type g =
      [=](const double& value1, const double& value2, OptionalJacobian<1,1> H1, OptionalJacobian<1, 1> H2) {
        if (H1)
          *H1 << value2;
        if (H2)
          *H2 << value1;
        return value1 * value2;
      };
  return double_(g, expression1, expression2);
}

double_ DivExpression(
    const double_& expression1, const double_& expression2) {
  typename double_::template BinaryFunction<double,double>::type g =
      [=](const double& value1, const double& value2, OptionalJacobian<1,1> H1, OptionalJacobian<1, 1> H2) {
        if (H1)
          *H1 << 1.0/value2;
        if (H2)
          *H2 << value1/(value2*value2);
        return value1 / value2;
      };
  return double_(g, expression1, expression2);
}

double_ gammaExpression(const double_& mu,
        const double_& Fz, const double_& U2);

double_ FyrExpression(const double_& mu,
    const double_& Fz, const double_& gamma,
    const double_& cAlpha, const double_& alpha);

double_ FyfExpression(const double_& mu,
    const double_& Fz,
    const double_& cAlpha, const double_& alpha);

double_ betadotExpression(const double_& mu,
    const double_& Fz, const double_& gamma,
    const double_& cAlpha, const double_& alpha,
    const double_& m, const double_& Vx,
    const double_& thetadd, const double_& Fyf,
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
