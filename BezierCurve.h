#pragma once
#ifndef BEZIERCURVE_H
#define BEZIERCURVE_H

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;


/**
 * @brief 三阶贝塞尔曲线类，在这个类中完成贝塞尔曲线的拟合
 */
class BezierCurve
{
public:
	BezierCurve(const std::vector<Vector2d>& _points);

	BezierCurve(const std::vector<Vector2d>& _points,const std::vector<Vector2d>& _tangent);
	~BezierCurve() {};
	std::vector<Eigen::Vector2d> outContralPoints()const;
private:
	void denseData(const std::vector<Vector2d>& _points);
	void doComputeTangent();
	void doComputeTangent_LSE();

	inline Eigen::Vector2d computeDerivative(const vector<Vector2d>& contralPoints,const double& t);
	inline Eigen::Vector2d computePrimePrime(const vector<Vector2d>& contralPoints, const double& t);
	inline Eigen::Vector2d pointAt(const vector<Vector2d>& controlPoints,const double& t);
	double simpson_3_8(const vector<Vector2d>& contralPoints, const double& leftTime, const double& rightTime);
	double adaptive_simpson_3_8(const vector<Vector2d>& contralPoints, const double& leftTime, const double& rightTime,const double& TOLERENCE = 1e-1);
	void doComputeLength(const vector<Vector2d>& contralPoints, const double& beginT=0.0, const double& endT=1.0);
	std::vector<Vector2d> bezierGenerate(vector<Vector2d>& points_, const int& begin, const int& end,vector<double>& parameter);
	std::vector<double> arcLengthParameterize(vector<Vector2d>& points_, const int& begin, const int& end);
	void reParameterize(const vector<Vector2d>& contralPoints, const vector<Vector2d>& points,const int& begin, const int& end, vector<double>& parameter);
	double computeMaxError(vector<Vector2d>& contraPoints, vector<Vector2d>& points_, const int& begin, const int& end,vector<double>& parameter,int& maxErrorIndex);
	vector<Vector2d> smoothData(vector<Vector2d>& points_, const int& begin, const int& end);
	void bezierCurveFitting(std::vector<Vector2d>& points_, const int& begin, const int& end,const double& TOLERENCE = 1e-2);
	bool bezierCurveFittingWithSmooth(std::vector<Vector2d>& points_, const int& begin, const int& end, const double& TOLERENCE = 5.0);

private:
	std::vector<Vector2d> points;
	std::vector<Vector2d> tangent;
	std::vector<Vector2d> ctrlPoints;
	//std::vector<double> arcLengthParameter;
	double m_length{-1};

};

#endif // !BEZIERCURVE_H

