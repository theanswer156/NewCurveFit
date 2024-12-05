#pragma once
#ifndef BSPLINE_H
#define BSPLINE_H
#include <iostream>
#include <Eigen/Core>
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace Eigen;
using namespace std;

class BSplineCurve {
public:
    BSplineCurve(const std::vector<Eigen::Vector2d>& controlPoints, const std::vector<double>& knotVector);

    void discreteCurve2Points();
    void BoehmKnotInsert();
    void OlsoKnotInsert();

    vector<double> outKnots() { return knotVector; }
    vector<Vector2d> outControlPoints() { return controlPoint; }
    vector<Vector2d> outCurvePoints() { return curvePoints; }
    vector<Vector2d> outPointsTangent() { return pointsTangent; }


    Vector2d PointAt(double& t);
    double basisFunction(const int& i, const double& t, const int& degree);
    Vector2d firstDiff(const double& t, const int& degree);
    Vector2d secondDiff(const double& t, const int& degree);
    void isKnotCorrect();
private:
    int degree = 0; //！B样条曲线的阶数而不是次数 次数等于阶数-1
    vector<Vector2d> controlPoint;
    vector<Vector2d> curvePoints;
    vector<Vector2d> pointsTangent;
    vector<double> knotVector;

};
#endif // BSPLINE_H


