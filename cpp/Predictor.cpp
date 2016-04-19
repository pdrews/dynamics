#include <gtsam/geometry/Pose2.h>
#include <functional>
#include <vector>
#include <cmath>
#include "Predictor.h"


/* x \in X == R^n
 * xdot \in Xdot == R^n
 * u \in U = R^p
 * theta \in Theta == R^m
 */

//n = 3
//p = 2
//m = 5

namespace gtsam{

//typedef Eigen::Matrix<double,4,100> MatrixBasis;
//typedef Eigen::Matrix<double,100,1> MatrixTheta;
//typedef Eigen::Matrix<double,4,1> MatrixX;
//typedef Eigen::Matrix<double,2,1> MatrixU;

/* Calculate basis functions for a specified X and U
 * X is expected to be:
 * roll
 * Ux
 * Uy
 * yaw_rate
 *
 * U is expected to be:
 * steering
 * throttle
 */
MatrixBasis Basis(MatrixX x, MatrixU u, double heading) {
  // For now, just some hard coded centroids
  // And for now, we will ignore the input u
  MatrixBasis b = MatrixBasis::Zero();
  Eigen::Matrix<double,THETA_DIM/X_DIM,1> basis;

  double alpha_f = 0;
  double alpha_r = 0;
  double steering = u(0);
  double throttle = u(1);
  double sd = sin(steering);

  if(x(1) > 0.1) {
    alpha_f = atan(x(2)/x(1) + 0.45*x(3)/x(1)) - steering;
    alpha_r = atan(x(2)/x(1) + 0.35*x(3)/x(1));
    basis(9) = x(2)/x(1)/40.0;
  }
  else {
    alpha_f = -steering;
    basis(9) = 0.0;
  }
  double kappa = (x(4)-(x(2)*(13.64/8.764)))/(x(2)*(13.64/8.764));

//  std::cout << "af " << alpha_f << " ar " << alpha_r << " steering " << steering << " throttle " << throttle << " sd " << sd << std::endl;

  basis(0) = throttle;
  basis(1) = x(1)/10.0;
  basis(2) = sd*tan(alpha_f)/1200.0;
  basis(3) = sd*tan(alpha_f) * std::abs(tan(alpha_f))/(1200.0*1200.0);
  basis(4) = sd*pow(tan(alpha_f),3)/(1200.0*1200.0*1200.0);
  basis(5) = x(3)*x(2)/25.0;
  basis(6) = x(3)/10.0;
  basis(7) = x(2)/10.0;
  basis(8) = sd;
  basis(10) = tan(alpha_f)/1400.0;
  basis(11) = tan(alpha_f) * std::abs(tan(alpha_f))/(1400.0*1400.0);
  basis(12) = pow(tan(alpha_f),3)/(pow(1400.0,3));
  basis(13) = tan(alpha_r)/(40.0);
  basis(14) = tan(alpha_r)*std::abs(tan(alpha_r))/(pow(40.0,2));
  basis(15) = pow(tan(alpha_r),3)/pow(40.0,3);
  basis(16) = x(3)*x(1)/50.0;
  basis(17) = x(0);
  basis(18) = x(0)*x(3);
  basis(19) = x(0)*x(1)/3.0;
  basis(20) = x(0)*x(1)*x(3)/5.0;
  basis(21) = heading/5.0;
  basis(22) = heading*x(0);
  basis(23) = pow(throttle,2);
  basis(24) = pow(throttle,3);
  basis(25) = kappa/(1.0+kappa);
  //basis(25) =
  //basis(25) =
  //basis(25) =

  b.block(0,0,1,THETA_DIM/X_DIM) = basis.transpose();
  b.block(1,THETA_DIM/X_DIM,1,THETA_DIM/X_DIM) = basis.transpose();
  b.block(2,2*THETA_DIM/X_DIM,1,THETA_DIM/X_DIM) = basis.transpose();
  b.block(3,3*THETA_DIM/X_DIM,1,THETA_DIM/X_DIM) = basis.transpose();
  b.block(4,4*THETA_DIM/X_DIM,1,THETA_DIM/X_DIM) = basis.transpose();

//  std::cout << "Basis Matrix: " << b << std::endl;
//  std::cout << "Small Matrix: " << basis << std::endl;

  return b;

}

typedef std::function<MatrixBasis(MatrixX, MatrixU, double)> BASIS;


  Predictor::Predictor(MatrixX x, MatrixU u, double heading, double dt) {
    x_ = x;
    BASIS basis(Basis);
    dt_basis_ = dt * basis(x, u, heading);
  }
  MatrixX Predictor::operator()(MatrixTheta theta, OptionalJacobian<X_DIM,THETA_DIM> H) {
    // calculate f(x,y) as sum of basis functions:
    // i.e., \sum w_i rbf(x,u)
    if (H) *H = dt_basis_;
    return x_ + dt_basis_*theta; // nxm * mx1
  }


}
