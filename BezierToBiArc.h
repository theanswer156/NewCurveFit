#pragma once
#ifndef BEZIER_TO_BIARC_H
#define BEZIER_TO_BIARC_H
#include <vector>
#include <Eigen/Core>
#include "BiArc.h"

using namespace Eigen;
using namespace std;


/*
* @brief 将Bezier曲线转换为BiArc曲线
* @param controlPoints Bezier曲线的控制点
* @param m_vecBiArc 存储BiArc曲线的容器
* @note 该类主要用于将Bezier曲线转换为BiArc曲线，并存储在m_vecBiArc中。
* @note 如果不在贝塞尔曲线曲率变化处将曲线进行二分，我们在进行拟合计算最大误差点的时候很可能会找到这个点，
*		与其后面在拟合的时候计算出来，不如先将曲率突变的点计算出来，并将贝塞尔曲线在此处分割，然后再做拟合。
***/
class BezierToBiArc
{
public:
	BezierToBiArc(const std::vector<Vector2d>& controlPoints);
	~BezierToBiArc() {};
	std::vector<BiArc> outBiArcs()const;
	inline Vector2d PointAt(const std::vector<Vector2d>& controlPoints, const double& t);
	inline Eigen::Vector2d computeDerivative(const vector<Vector2d>& contralPoints, const double& t);
	Eigen::Vector2d innerPoint(const Vector2d& p1, const Vector2d& p2, const Vector2d& p3);
	Eigen::Vector2d lineIntersection(const Vector2d& p1, const Vector2d& t1, const Vector2d& p2, const Vector2d& t2);
	Eigen::Vector2d computeCircleCenter(const Vector2d& p1, const Vector2d& t, const Vector2d& p2);
	void splitBezierByCurvature(const std::vector<Vector2d>& controlPoints);
	void fromCurvesToBiArc();
	void singleBezierToBiArc(const std::vector<Vector2d>& controlPoints, const double& beginTime = 0.0, const double& endTime = 1.0,const double& TOLERENCE = 2.0);
	void fromBezierToBiArc(const std::vector<Vector2d>& controlPoints, const double& beginTime = 0.0, const double& endTime = 1.0,const double& TOLERENCE = 2.0);
	void fromBezierToBiArc_IP(const std::vector<Vector2d>& controlPoints, const double& beginTime = 0.0, const double& endTime = 1.0, const double& TOLERENCE = 2.0);
	inline double vectorAngle(const Vector2d& v1, const Vector2d& v2);
	inline double crossProduct(const Vector2d& v1, const Vector2d& v2);
	vector<double> computeMaxError(BiArc& biarc, const vector<Vector2d>& controlPoints, const double& beginTime, const double& endTime, double& maxErrorIndex);
	double getCurvateMutations(const std::vector<Vector2d>& controlPoints);
private:
	std::vector<BiArc> m_vecBiArc = {};
	std::vector<std::vector<Vector2d>> m_vecControlPoints = {};
	double m_TOLERENCE = 2.0;
};

#endif // !BEZIER_TO_BIARC_H