/*
 * testFrame.cpp
 *
 *  Created on: Feb 13, 2014
 *      Author: cbeall3
 */

#include <vo/Frame.h>
#include <CppUnitLite/TestHarness.h>
#include <gtsam/base/Testable.h>
#include <gtsam/nonlinear/expressions.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/ExpressionFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/geometry/Point3.h>
#include "../PredictorExpression.h"


using namespace std;
using namespace gtsam;

typedef Expression<Point3> Point3_;
typedef Expression<Pose3> Pose3_;
typedef Expression<Rot3> Rot3_;

/* ************************************************************************* */
// Many Leaves
TEST(Vec1, test) {
  Rot3_ R1(1), R2(2);
  Rot3_ R3 = R1 * R2;

  const Key key1(67), key2(68);
  const Point3_ sum_ = Point3_(key1) + Point3_(Point3(1, 1, 1));

  Expression<double> e1 = Expression<double>(key1) + Expression<double>(key2);
  Expression<double> e2 = Expression<double>(key1) - Expression<double>(key2);
  Expression<double> e3 = MulExpression(Expression<double>(key1), Expression<double>(key2));
  Expression<double> e4 = DivExpression(Expression<double>(key1), Expression<double>(key2));

  Values values;

//  Expression<Point3> a1_ = Expression<Point3>(1);
//  Expression<Point3> a2_ = Expression<Point3>(2);
//  Expression<Point3> a3_ = a1_+a2_;
//  Expression<Point2> a4_ = a1_ - a2_;

//  Vector1 a1,a2;
//  a1 << 1;
//  a2 << 2;
//  Vector2 expected(1,2);
//  Expression<Vector1> e1(10);
//  Expression<Vector1> e2(11);
//  values.insert(10, a1);
//  values.insert(11, a2);
//  Expression<Vector2> exp = StackExpression2<Vector2,Vector1>(e1, e2);
//
//  cout << " Values 10 " << values.at<Vector1>(10) << endl;
//  cout << " Values 11 " << values.at<Vector1>(11) << endl;
//  cout << " Vector a1 " << a1 << endl;
//  cout << " Vector a2 " << a2 << endl;
//
//  EXPECT(assert_equal(expected, exp.value(values)));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
