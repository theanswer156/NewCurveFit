#pragma once
#ifndef CUBICBEZIER_H
#define CUBICBEZIER_H

#include "BezierCurve.h"
#include <memory>

class CubicBezier:public BezierCurve
{
public:
	CubicBezier() = default;
	CubicBezier(const std::vector<Vector2d>& controlPoints);
	CubicBezier(const std::vector<Vector2d>& contralPoints, const std::pair<double, double>& range);

	const vector<CubicBezier> splitDXF(const std::string& filename);
	vector<CubicBezier>& transSplineToBezier(std::ifstream& infile);

	inline friend void operator<<(std::ostream& os, const CubicBezier& cb);

	~CubicBezier() = default;

private:
	std::pair<double, double> m_vecRange{-1,-1};
};

#endif // !CUBICBEZIER_H


