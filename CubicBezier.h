#pragma once
#ifndef CUBICBEZIER_H
#define CUBICBEZIER_H

#include "BezierCurve.h"
//#include <memory>
#include <ranges>

/*
*@brief 三次贝塞尔曲线类，继承自BezierCurve类
* @param controlPoints 控制点集合
* @param range 范围(这里不再默认是[0,1],而是由任意值,由此,在绘制这个类的图形时，要根据范围进行变化)
**/
class CubicBezier:public BezierCurve
{
public:
	CubicBezier() = default;
	CubicBezier(const std::vector<Vector2d>& controlPoints);
	CubicBezier(const std::vector<Vector2d>& contralPoints, const std::pair<double, double>& range);
	CubicBezier(const std::ranges::subrange<vector<Vector2d>::const_iterator> controlPoints, const std::pair<double, double>& range);

	const vector<CubicBezier> splitDXF(const std::string& filename);
	vector<CubicBezier> transSplineToBezier(std::ifstream& infile);
	bool isUniform(const std::vector<double>& knotVec);
	void uniformSplineToBezier(const std::vector<double>& knotVec,const std::vector<Vector2d>& controlPoints, std::vector<CubicBezier>& bezierVec);
	void NonuniformSplineToBezier(const std::vector<double>& knotVec, const std::vector<Vector2d>& controlPoints, std::vector<CubicBezier>& bezierVec);
	Eigen::Vector2d DeRecursion(const vector<double>& knotVec, const vector<Vector2d>& controlPoints, const int& i, const int& r, const int& t);
	inline friend void operator<<(std::ostream& os, const CubicBezier& cb);

	~CubicBezier() = default;

private:
	std::pair<double, double> m_vecRange{-1,-1};
};

#endif // !CUBICBEZIER_H


