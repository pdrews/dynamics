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

double_ gammaExpression(const double_& mu,
        const double_& Fz, const double_& U2){
  double_ muFz = MulExpression(mu,Fz);
  double_ e1 = MulExpression(muFz, muFz) - MulExpression(U2,U2);
  return DivExpression(e1,muFz);
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
  double_ e1 = tanhExpression(MulExpression(cAlpha,alpha));
  return(-1*MulExpression(gamma,MulExpression(mu,MulExpression(Fz,e1))));
//  return(-1*gamma*mu*Fz*tanhExpression(cAlpha*alpha));
}

double_ FyfExpression(const double_& mu,
    const double_& Fz,
    const double_& cAlpha, const double_& alpha) {
  double_ e1 = tanhExpression(MulExpression(cAlpha,alpha));
  return(-1*MulExpression(mu,MulExpression(Fz,e1)));
//  return(-1*mu*Fz*tanhExpression(cAlpha*alpha));
}

double_ betadotExpression(const double_& mu,
    const double_& Fz, const double_& gamma,
    const double_& cAlpha, const double_& alpha,
    const double_& m, const double_& Vx,
    const double_& thetadd, const double_& Fyf,
    const double_& Fyr) {
  return DivExpression((Fyf + Fyr),MulExpression(m,Vx)) - thetadd;
}

double_ VxdExpression(const double_& U2,
    const double_& Fyf, const double_& U1,
    const double_& m, const double_& Vx,
    const double_& beta, const double_& thetad) {
  double_ e1 = MulExpression(Vx,MulExpression(beta,thetad));
  return DivExpression((U2 - MulExpression(Fyf,sinExpression(U1))),m)+ e1;
//  return ((U2 - (Fyf*sinExpression(U1)))/m)+ (Vx*beta*thetad);
}

double_ thetaddExpression(const double_& a,
    const double_& Fyf, const double_& b,
    const double_& Fyr, const double_& Iz) {
  return DivExpression(MulExpression(a,Fyf) - MulExpression(b,Fyr),Iz);
//  return ((a*Fyf) - (b*Fyr))/Iz;
}

double_ alphafExpression(const double_& beta,
    const double_& a, const double_& Vx,
    const double_& thetad, const double_& delta) {
  return arctanExpression(beta + MulExpression(DivExpression(a,Vx), thetad)) - delta;
}

double_ alpharExpression(const double_& beta,
    const double_& a, const double_& Vx,
    const double_& thetad) {
  return arctanExpression(beta - MulExpression(DivExpression(a,Vx), thetad));
//  return arctanExpression(beta - ((a/Vx) * thetad));
}
}
