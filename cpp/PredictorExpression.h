/**
 * @file Predictor.h
 * @brief Dynamics Predictor
 * @author Paul Drews
 */

#pragma once

#include <gtsam_unstable/base/dllexport.h>
#include <gtsam/geometry/Pose3.h>
#include <math.h>

namespace gtsam {

/**
 * Create an expression out of a linear function f:T->A with (constant) Jacobian dTdA
 * Specialize for stacking scalars to vector
 * T is the return vector type, A is the scalar type of the expression
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
Expression<double> cosExpression(const Expression<double>& expression) {
  // Use lambda to calculate the Jacobian
  typename Expression<double>::template UnaryFunction<double>::type g =
      [=](const double& value, OptionalJacobian<1,1> H) {
        if (H)
          *H << -sin(value);
        return cos(value);
      };
  return Expression<double>(g, expression);
}

/**
 * Create a cosine factor
 */
Expression<double> sinExpression(const Expression<double>& expression) {
  // Use lambda to calculate the Jacobian
  typename Expression<double>::template UnaryFunction<double>::type g =
      [=](const double& value, OptionalJacobian<1,1> H) {
        if (H)
          *H << cos(value);
        return sin(value);
      };
  return Expression<double>(g, expression);
}

Expression<double> arctanExpression(const Expression<double>& expression) {
  // Use lambda to calculate the Jacobian
  typename Expression<double>::template UnaryFunction<double>::type g =
      [=](const double& value, OptionalJacobian<1,1> H) {
        if (H)
          *H << 1/((value*value)+1);
        return atan(value);
      };
  return Expression<double>(g, expression);
}

Expression<double> sqrtExpression(const Expression<double>& expression) {
  // Use lambda to calculate the Jacobian
  typename Expression<double>::template UnaryFunction<double>::type g =
      [=](const double& value, OptionalJacobian<1,1> H) {
        if (H) {
          assert(value > 0);
          *H << 1/(2*sqrt(value));
        }
        return sqrt(value);
      };
  return Expression<double>(g, expression);
}

} // \namespace gtsam
