#include "dynamics.h"
#include "PredictorExpression.h"
//#include "Predictor.h"

#include <gtsam/nonlinear/Expression.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/ExpressionFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/PriorFactor.h>


#include <boost/tokenizer.hpp>
#include "boost/program_options.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>

using namespace gtsam;
using namespace std;
using symbol_shorthand::T;
using symbol_shorthand::V;
using symbol_shorthand::B;
using symbol_shorthand::P;
int main()
{
  namespace po = boost::program_options;
  //Load options from config file
  double lambda, muv, massv, av, bv, Izv, muSigmav, massSigmav, aSigmav, bSigmav, IzSigmav, c1v, c2v, c3v, calphav;
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produce help message")
      ("mu", po::value<double>(&muv), "")
      ("mass", po::value<double>(&massv), "")
      ("a", po::value<double>(&av), "")
      ("b", po::value<double>(&bv), "")
      ("Iz", po::value<double>(&Izv), "")
      ("c1", po::value<double>(&c1v), "")
      ("c2", po::value<double>(&c2v), "")
      ("c3", po::value<double>(&c3v), "")
      ("calpha", po::value<double>(&calphav), "")
      ("lambda", po::value<double>(&lambda), "")
      ("muSigma", po::value<double>(&muSigmav), "")
      ("massSigma", po::value<double>(&massSigmav), "")
      ("aSigma", po::value<double>(&aSigmav), "")
      ("bSigma", po::value<double>(&bSigmav), "")
      ("IzSigma", po::value<double>(&IzSigmav), "")
  ;

  po::variables_map options;
  std::ifstream file("../conf.txt");
  po::store(po::parse_config_file(file, desc), options);
  file.close();
  po::notify(options);

  cout << muv << endl << massv << endl << av << endl <<
      bv << endl << Izv << endl << muSigmav << endl <<
      massSigmav << endl << aSigmav << endl <<
      bSigmav << endl << IzSigmav << endl;

  //Load up all of our data
  ifstream data;
  data.open("../data.txt");
  string line;
  int numSamples = 0;
  while(getline(data, line)) {numSamples++;}
  cout << "Found " << numSamples << " Samples" << endl;
  Eigen::Matrix<double,X_DIM,Eigen::Dynamic> X(X_DIM,numSamples);
  Eigen::Matrix<double,U_DIM,Eigen::Dynamic> U(U_DIM,numSamples);
  Eigen::Matrix<double,1,Eigen::Dynamic> heading(1,numSamples);
  Eigen::Matrix<double,1,Eigen::Dynamic> betam(1,numSamples);
  Eigen::Matrix<double,1,Eigen::Dynamic> Vxm(1,numSamples);
  Eigen::Matrix<double,1,Eigen::Dynamic> ThetaDm(1,numSamples);

  data.close();
  data.open("../data.txt");
  int sampleNum = 0;
  while(getline(data, line)) {
    double left, right;
    boost::char_separator<char> sep(",");
    boost::tokenizer<boost::char_separator<char>> toks(line,sep);
    boost::tokenizer<boost::char_separator<char>>::iterator tok = toks.begin();
    //Roll
    X(0,sampleNum) = atof(tok->c_str());
    tok++;
    //Ux
    X(1,sampleNum) = atof(tok->c_str());
    tok++;
    //Uy
    X(2,sampleNum) = atof(tok->c_str());
    tok++;
    //Yaw rate
    X(3,sampleNum) = atof(tok->c_str());
    tok++;
    //Yaw
    heading(0,sampleNum) = atof(tok->c_str());
    tok++;
    //Throttle
    U(0,sampleNum) = atof(tok->c_str());
    tok++;
    //Steering
    U(1,sampleNum) = atof(tok->c_str());
    left = atof(tok->c_str());
    tok++;
    right = atof(tok->c_str());
    X(4,sampleNum) = (left+right) / 2.0;

    betam(0,sampleNum) = atan2(X(2,sampleNum),X(1,sampleNum));
    ThetaDm(0,sampleNum) = X(3,sampleNum);
    Vxm(0,sampleNum) = X(1,sampleNum);

    sampleNum++;
  }
  cout << "Finished reading " << sampleNum << " samples" << endl;
  double dt = 1.0/200.0;

  //Define expressions for all of our parameters
  double_ Calpha('p',0);    //About 10000???
  double_ mu('p',1);        //About 0.5??
  double_ mass('p',2);      //About 20kg
  double_ a('p',3);         //About .45m
  double_ b('p',4);         //About .35m
  double_ Iz('p', 5);       //About 2??????
  //Parameters for fitting steering and throttle
  double_ c1('p',6); //throttle
  double_ c2('p',7); //steering
  double_ c3('p',8); //throttle offset

  //Magic tire model
  //Derived parameters
  double_ Fzf = (9.8 * 0.5) * Mul_(Div_(b,(a+b)), mass);
  double_ Fzr = (9.8 * 0.5) * Mul_(Div_(a,(a+b)), mass);
//  double_ Fzf(10.0);
//  double_ Fzr(350.0);

  NonlinearFactorGraph graph;
  Values vals;

//  PriorFactor<double> priorCalpha(P(0), 0.5, noiseModel::Diagonal::Sigmas(Matrix1::Constant(0.3)));
  PriorFactor<double> priormu(P(1), muv, noiseModel::Diagonal::Sigmas(Matrix1::Constant(muSigmav)));
  PriorFactor<double> priormass(P(2), massv, noiseModel::Diagonal::Sigmas(Matrix1::Constant(massSigmav)));
  PriorFactor<double> priora(P(3), av, noiseModel::Diagonal::Sigmas(Matrix1::Constant(aSigmav)));
  PriorFactor<double> priorb(P(4), bv, noiseModel::Diagonal::Sigmas(Matrix1::Constant(bSigmav)));
  PriorFactor<double> priorIz(P(5), Izv, noiseModel::Diagonal::Sigmas(Matrix1::Constant(IzSigmav)));
//  PriorFactor<double> priorc1(P(6), 0.35, noiseModel::Diagonal::Sigmas(Matrix1::Constant(0.1)));
//  PriorFactor<double> priorc2(P(7), 0.5, noiseModel::Diagonal::Sigmas(Matrix1::Constant(0.3)));
//  PriorFactor<double> priorc3(P(8), 0.5, noiseModel::Diagonal::Sigmas(Matrix1::Constant(0.3)));

//  graph.add(priorCalpha);
  graph.add(priormu);
  graph.add(priormass);
  graph.add(priora);
  graph.add(priorb);
  graph.add(priorIz);
//  graph.add(priorc1);
//  graph.add(priorc2);
//  graph.add(priorc3);

  //Put in initial guesses
  vals.insert(P(0), calphav);
  vals.insert(P(1), muv);
  vals.insert(P(2), massv);
  vals.insert(P(3), av);
  vals.insert(P(4), bv);
  vals.insert(P(5), Izv);

  vals.insert(P(6), c1v);
  vals.insert(P(7), c2v);
  vals.insert(P(8), c3v);


  //Iterate through each pair in our data, putting it in our FG
  Matrix1 sigma1 = Matrix1::Constant(0.1);
  Matrix1 sigma2 = Matrix1::Constant(0.01);

  noiseModel::Diagonal::shared_ptr R = noiseModel::Diagonal::Sigmas(sigma1);
  noiseModel::Diagonal::shared_ptr R2 = noiseModel::Diagonal::Sigmas(sigma2);


  for(int i=0; i<numSamples-2; i++) {
//  for(int i=0; i<2; i++) {
    //put vx, beta, thetadot in the graph
//    vals.insert(B(i), betam(0,i));
//    vals.insert(T(i), ThetaDm(0,i));
//    vals.insert(V(i), Vxm(0,i));
    //Prior factors
//    PriorFactor<double> priorB(B(i), betam(0,i), R);
//    PriorFactor<double> priorT(T(i), ThetaDm(0,i), R);
//    PriorFactor<double> priorV(V(i), Vxm(0,i), R);
//    graph.add(priorB);
//    graph.add(priorT);
//    graph.add(priorV);
    double_ beta_(betam(0,i)), thetad_(ThetaDm(0,i)), vx_(Vxm(0,i));
//    double_ beta_('b',i), thetad_('t',i), vx_('v',i);
    double_ steering_ = U(1,i)*c2;
      double_ throttle_ = Mul_(mu,Mul_(Fzr,tanh_(Mul_(U(0,i),c1)))) - Vxm(0,i)*c3;
//        double_ throttle_ = U(0,i)*c1;// + Vxm(0,i)*c3;

    if(std::abs(X(1,i) - X(1,i+1)) > 15) {
      cout << X(1,i) << " " << X(1,i+1) << endl;
      cout << "Skipping a sample between bag files" << endl;
      continue;
    }
    //Add prediction factor requirements
    double_ gamma = gammaExpression(mu, Fzr, throttle_);
    double_ alphaf = alphafExpression(beta_, a, vx_, thetad_, steering_);
    double_ alphar = alpharExpression(beta_, a, vx_, thetad_);
    double_ Fyr = FyrExpression(mu, Fzr, gamma, Calpha, alphar);
    double_ Fyf = FyfExpression(mu, Fzr, Calpha, alphaf);

    //Add prediction factors
    double_ betadot = betadotExpression(mu, mass, vx_, thetad_, Fyf, Fyr);
    double_ vxdot = VxdExpression(steering_, Fyf, throttle_, mass, vx_, beta_, thetad_);
    double_ thetaddot = thetaddExpression(a, Fyf, b, Fyr, Iz);

    graph.addExpressionFactor<double>(R, betam(0,i+1) - betam(0,i), (dt*betadot));
    graph.addExpressionFactor<double>(R, Vxm(0,i+1) - Vxm(0,i), (dt*vxdot));
    graph.addExpressionFactor<double>(R2, ThetaDm(0,i+1) - ThetaDm(0,i), (dt*thetaddot));

//    cout << "Expression for Beta " << i << endl;
//    betadot.print("Beta: ");

  }
//  graph.print("Our factor graph:\n");
//  vals.print("InitialValues\n");


  //Solve the FG and spit out the thetas
//  MatrixTheta Theta = MatrixTheta::Zero();
//  vals.insert(50, Theta);
//  //  graph.print("Factor Graph:\n");
  LevenbergMarquardtParams params;
  params.lambdaInitial = lambda;
//  params.verbosity = LevenbergMarquardtParams::Verbosity::DELTA;
  params.verbosityLM = LevenbergMarquardtParams::VerbosityLM::SUMMARY;
  LevenbergMarquardtOptimizer optimizer(graph, vals, params);
  optimizer.iterate();



//  GaussNewtonParams params;
//  params.verbosity = GaussNewtonParams::Verbosity::DELTA;
//  GaussNewtonOptimizer optimizer(graph, vals, params);

  Values result = optimizer.optimize();
  cout << "We did " << optimizer.iterations() << " Iterations" << endl;
//  optimizer.print("Optimizer:");
  cout << "mu: " << result.at<double>(P(1)) << endl;
  cout << "mass: " << result.at<double>(P(2)) << endl;
  cout << "a: " << result.at<double>(P(3)) << endl;
  cout << "b: " << result.at<double>(P(4)) << endl;
  cout << "Calpha: " << result.at<double>(P(0)) << endl;
  cout << "Iz: " << result.at<double>(P(5)) << endl;
  cout << "c1: " << result.at<double>(P(6)) << endl;
  cout << "c2: " << result.at<double>(P(7)) << endl;
//  cout << "c3: " << result.at<double>(P(8)) << endl;

      result.print("Final Result:\n");
//
//  // Get the average and mean squared error
//  Eigen::Matrix<double,X_DIM,Eigen::Dynamic> der_pred(X_DIM,numSamples);
//  Eigen::Matrix<double,X_DIM,Eigen::Dynamic> der_actual(X_DIM,numSamples);
//
    Eigen::Matrix<double,3,1> squared_error = Eigen::Matrix<double,3,1>::Zero();
//  MatrixX squared_error_der = MatrixX::Zero();
//  MatrixBasis J;
//  MatrixTheta thetaEst = result.at<MatrixTheta>(50);
  ofstream sed, sedp, cov;
  sed.open("../sed.txt");
  sedp.open("../sedp.txt");
  int skips = 0;
  for(int i=0; i<numSamples-1; i++) {
    if(std::abs(X(1,i) - X(1,i+1)) > 15) {
      cout << X(1,i) << " " << X(1,i+1) << endl;
      cout << "Skipping a sample between bag files" << endl;
      continue;
    }
    //Setup prediction
    double_ beta_(betam(0,i)), thetad_(ThetaDm(0,i)), vx_(Vxm(0,i));
    double_ steering_ = U(1,i)*c2;
    double_ throttle_ = U(0,i)*c1;
    //Add prediction factor requirements
    double_ gamma = gammaExpression(mu, Fzr, throttle_);
    double_ alphaf = alphafExpression(beta_, a, vx_, thetad_, steering_);
    double_ alphar = alpharExpression(beta_, a, vx_, thetad_);
    double_ Fyr = FyrExpression(mu, Fzr, gamma, Calpha, alphar);
    double_ Fyf = FyfExpression(mu, Fzr, Calpha, alphaf);
    //Add prediction factors
    double_ betadot = betadotExpression(mu, mass, vx_, thetad_, Fyf, Fyr);
    double_ vxdot = VxdExpression(steering_, Fyf, throttle_, mass, vx_, beta_, thetad_);
    double_ thetaddot = thetaddExpression(a, Fyf, b, Fyr, Iz);

    double betapred = betam(0,i) + dt*betadot.value(result);
    double vxpred = Vxm(0,i) + dt*vxdot.value(result);
    double thetaddotpred = ThetaDm(0,i) + dt*thetaddot.value(result);

    squared_error(0) += (betapred-betam(0,i+1))*(betapred-betam(0,i+1));
    squared_error(1) += (vxpred-Vxm(0,i+1))*(vxpred-Vxm(0,i+1));
    squared_error(2) += (thetaddotpred-ThetaDm(0,i+1))*(thetaddotpred-ThetaDm(0,i+1));

  }
  squared_error = squared_error.array() / ((numSamples-skips)-1);
  cout << "Squared error is \n" << squared_error << endl;

}
