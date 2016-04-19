/*
 * testFrame.cpp
 *
 *  Created on: Feb 13, 2014
 *      Author: cbeall3
 */

#include <vo/Frame.h>
#include <CppUnitLite/TestHarness.h>
#include <gtsam/base/Testable.h>
#include <gtsam/nonlinear/Expression.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/ExpressionFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include "../PredictorExpression.h"


using namespace std;
using namespace gtsam;

/* ************************************************************************* */
// Many Leaves
TEST(Stack, Exp2) {
  Values values;
  Vector1 a1,a2;
  a1 << 1;
  a2 << 2;
  Vector2 expected(1,2);
  Expression<Vector1> e1(10);
  Expression<Vector1> e2(11);
  values.insert(10, a1);
  values.insert(11, a2);
  Expression<Vector2> exp = StackExpression2<Vector2,Vector1>(e1, e2);

  cout << " Values 10 " << values.at<Vector1>(10) << endl;
  cout << " Values 11 " << values.at<Vector1>(11) << endl;
  cout << " Vector a1 " << a1 << endl;
  cout << " Vector a2 " << a2 << endl;

  EXPECT(assert_equal(expected, exp.value(values)));
}

/* ************************************************************************* */
// Many Leaves
TEST(Stack, Exp2double) {
  Values values;
  double a1(1),a2(2);
  Vector2 expected(1,2);
  Expression<double> e1((Key)10);
  Expression<double> e2((Key)11);
  values.insert(10, a1);
  values.insert(11, a2);
  Expression<Vector2> exp = StackExpression2<Vector2,double>(e1, e2);

  cout << " Values 10 " << values.at<double>(10) << endl;
  cout << " Values 11 " << values.at<double>(11) << endl;
  cout << " Vector a1 " << a1 << endl;
  cout << " Vector a2 " << a2 << endl;

  EXPECT(assert_equal(expected, exp.value(values)));
}

TEST(Stack, Exp3) {
  Values values;
  Vector1 a1,a2,a3;
  a1 << 1;
  a2 << 2;
  a3 << 3;
  Vector3 expected(1,2,3);
  Expression<Vector1> e1(10);
  Expression<Vector1> e2(11);
  Expression<Vector1> e3(12);
  values.insert(10, a1);
  values.insert(11, a2);
  values.insert(12, a3);
  Expression<Vector3> exp = StackExpression3<Vector3,Vector1>(e1, e2, e3);

  EXPECT(assert_equal(expected, exp.value(values)));
}

TEST(Stack, cosExpression) {
  Values values;
  double a(1),b(0),c(2);

  Expression<double> e1((Key)10);
  Expression<double> e2((Key)11);
  Expression<double> e3((Key)12);
  values.insert(10, a);
  values.insert(11, b);
  values.insert(12, c);
  Expression<double> exp = cosExpression(e1);
  Matrix1 Hexp;
  std::vector<Matrix> H(1);
  Hexp << -sin(1);
  EXPECT(assert_equal(cos(1), exp.value(values, H)));
  EXPECT(assert_equal(Hexp, H[0]))

  Hexp << -sin(0);
  Expression<double> exp2 = cosExpression(e2);
  EXPECT(assert_equal(cos(0), exp2.value(values, H)));
  EXPECT(assert_equal(Hexp, H[0]))

  Hexp << -sin(2);
  Expression<double> exp3 = cosExpression(e3);
  EXPECT(assert_equal(cos(2), exp3.value(values, H)));
  EXPECT(assert_equal(Hexp, H[0]))

}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
