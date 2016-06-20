#include <gtsam/geometry/Pose2.h>
#include <functional>
#include <vector>
#include <cmath>
#include <gtsam/nonlinear/Expression.h>
#include "PredictorExpression.h"


/* x \in X == R^n
 * xdot \in Xdot == R^n
 * u \in U = R^p
 * theta \in Theta == R^m
 */

//n = 3
//p = 2
//m = 5

/* Parameters
 * Calpha
 * Fz
 * mu
 * mass
 * a
 * b
 */

/* States
 * derivative of sideslip
 */

namespace gtsam{

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

double_ tanh_(const double_& expression) {
  // Use lambda to calculate the Jacobian
  typename double_::template UnaryFunction<double>::type g =
      [=](const double& value, OptionalJacobian<1,1> H) {
        if (H)
          *H << pow(1/(cosh(value)),2);
        return tanh(value);
      };
  return double_(g, expression);
}

double_ sqrt_(const double_& expression) {
  // Use lambda to calculate the Jacobian
  typename double_::template UnaryFunction<double>::type g =
      [=](const double& value, OptionalJacobian<1,1> H) {
        if (H) {
          assert(value > 0);
          //if (value < 0) {
          //  std::cout << "Val too high!! " << value << std::endl;
         // }
          *H << 1/(2*sqrt(value));
        }
        return sqrt(value);
      };
  return double_(g, expression);
}

double_ Mul_(
    const double_& expression1, const double_& expression2) {
  typename double_::template BinaryFunction<double,double>::type g =
      [=](const double& x, const double& y, OptionalJacobian<1,1> H1, OptionalJacobian<1, 1> H2) {
        if (H1)
          *H1 << y;
        if (H2)
          *H2 << x;
        return x * y;
      };
  return double_(g, expression1, expression2);
}

double_ Div_(
    const double_& expression1, const double_& expression2) {
  typename double_::template BinaryFunction<double,double>::type g =
      [=](const double& x, const double& y, OptionalJacobian<1,1> H1, OptionalJacobian<1, 1> H2) {
        if (H1)
          *H1 << 1.0/y;
        if (H2)
          *H2 << x/(y*y);
        return x / y;
      };
  return double_(g, expression1, expression2);
}

double_ gammaExpression(const double_& mu,
        const double_& Fz, const double_& U2){
  double_ muFz = Mul_(mu,Fz);
  double_ e1 = Mul_(muFz, muFz) - Mul_(U2,U2);
  return Div_(sqrt_(e1),muFz);
//  return (((mu*Fz) * (mu*Fz)) - (U2*U2))/(mu*Fz);
}

//double_ alphaslExpression(const double_& mu,
//    const double_& Fz, const double_& gamma,
//    const double_& cAlpha){
//  return(arctanExpression((3*gamma*mu*Fz)/cAlpha));
//}

double_ FyrExpression(const double_& mu,
    const double_& Fz, const double_& gamma,
    const double_& cAlpha, const double_& alpha) {
  double_ e1 = tanh_(Mul_(cAlpha,alpha));
  return(-1*Mul_(gamma,Mul_(mu,Mul_(Fz,Div_(e1,gamma)))));
//  return(-1*gamma*mu*Fz*tanhExpression(cAlpha*alpha));
}

double_ FyfExpression(const double_& mu,
    const double_& Fz,
    const double_& cAlpha, const double_& alpha) {
  double_ e1 = tanh_(Mul_(cAlpha,alpha));
  return(-1.0*Mul_(mu,Mul_(Fz,e1)));
//  return(-1*mu*Fz*tanhExpression(cAlpha*alpha));
}

double_ betadotExpression(const double_& mu,
    const double_& m, const double_& Vx,
    const double_& thetad, const double_& Fyf,
    const double_& Fyr) {
  return Div_((Fyf + Fyr),Mul_(m,Vx)) - thetad;
}

double_ VxdExpression(const double_& steering,
    const double_& Fyf, const double_& throttle,
    const double_& m, const double_& Vx,
    const double_& beta, const double_& thetad) {
  double_ e1 = Mul_(Vx,Mul_(beta,thetad));
  return Div_((steering - Mul_(Fyf,sinExpression(throttle))),m)+ e1;
//  return ((U2 - (Fyf*sinExpression(U1)))/m)+ (Vx*beta*thetad);
}

double_ thetaddExpression(const double_& a,
    const double_& Fyf, const double_& b,
    const double_& Fyr, const double_& Iz) {
  return Div_(Mul_(a,Fyf) - Mul_(b,Fyr),Iz);
//  return ((a*Fyf) - (b*Fyr))/Iz;
}

double_ alphafExpression(const double_& beta,
    const double_& a, const double_& Vx,
    const double_& thetad, const double_& delta) {
  return arctanExpression(beta + Mul_(Div_(a,Vx), thetad)) - delta;
}

double_ alpharExpression(const double_& beta,
    const double_& a, const double_& Vx,
    const double_& thetad) {
  return arctanExpression(beta - Mul_(Div_(a,Vx), thetad));
//  return arctanExpression(beta - ((a/Vx) * thetad));
}
}
