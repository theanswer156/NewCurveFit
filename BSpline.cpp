#include "BSpline.h"



BSplineCurve::BSplineCurve(const std::vector<Eigen::Vector2d>& controlPoints, const std::vector<double>& knotVector)
	: controlPoint(controlPoints), knotVector(knotVector), degree(knotVector.size() - controlPoints.size() - 1)
{
	isKnotCorrect();
	discreteCurve2Points();
}

/*
*@brief 将B样条曲线离散化为一系列点
*
***/
void BSplineCurve::discreteCurve2Points()
{
	int n = controlPoint.size();
	double startT = *(knotVector.begin()+degree-1);
	double endT = *(knotVector.begin() + n);

	double change_range = endT - startT;
	//! 通过设置步长来控制离散化精度
	//! 最小取1/0.005(200个点)，大于200则取2e-2
	double step = change_range > 4.0 ? 5e-3 : 1e-3;

	double t = startT;
	while (t <= endT)
	{
		Eigen::Vector2d point = PointAt(t);
		Eigen::Vector2d tangent = firstDiff(t, degree).normalized();

		curvePoints.emplace_back(point);
		pointsTangent.emplace_back(tangent);
		t += step;
	}
}

Vector2d BSplineCurve::PointAt(double& t)
{
	Eigen::Vector2d ans(0.0, 0.0);
	for (int i = 0; i < controlPoint.size(); i++)
	{
		double basis = basisFunction(i, t, degree);
		ans += basis * controlPoint[i];
	}
	return ans;
}

double BSplineCurve::basisFunction(const int& i, const double& t, const int& _degree)
{
	if (!_degree)
	{
		return t >= knotVector[i] && t < knotVector[i + 1] ? 1.0 : 0.0;
	}

	double denominator1 = knotVector[i + _degree] - knotVector[i];
	double denominator2 = knotVector[i + _degree + 1] - knotVector[i + 1];

	double part1 = denominator1 > 0 ? (t - knotVector[i]) * basisFunction(i, t, _degree - 1) / denominator1 : 0.0;
	double part2 = denominator2 > 0 ? (knotVector[i + _degree + 1] - t) * basisFunction(i + 1, t, _degree - 1) / denominator2 : 0.0;

	return part1 + part2;
}

Vector2d BSplineCurve::firstDiff(const double& t, const int& _degree)
{
	Eigen::Vector2d ans(0.0, 0.0);

	for (int j = 0; j < controlPoint.size(); j++)
	{
		double basisDiff = 0.0;
		if (_degree > 0)
		{
			double denominator1 = knotVector[j + _degree] - knotVector[j];
			double denominator2 = knotVector[j + _degree + 1] - knotVector[j + 1];
			basisDiff += denominator1 > 0 ? _degree * basisFunction(j, t, _degree - 1) / denominator1 : 0.0;
			basisDiff += denominator2 > 0 ? -_degree * basisFunction(j + 1, t, _degree - 1) / denominator2 : 0.0;
		}
		ans += basisDiff * controlPoint[j];
	}

	return ans;
}

Vector2d BSplineCurve::secondDiff(const double& t, const int& _degree)
{
	Eigen::Vector2d ans(0.0, 0.0);

	for (int j = 0; j < controlPoint.size(); j++)
	{
		double basisSecondDiff = 0.0;
		if (_degree > 1)
		{
			double denominator1_1 = knotVector[j + _degree] - knotVector[j];
			double denominator1_2 = knotVector[j + _degree - 1] - knotVector[j];
			double denominator2_1 = knotVector[j + _degree + 1] - knotVector[j + 1];
			double denominator2_2 = knotVector[j + _degree] - knotVector[j + 1];

			if (denominator1_1 > 0 && denominator1_2 > 0)
			{
				basisSecondDiff += _degree * (_degree - 1) * basisFunction(j, t, _degree - 2) / (denominator1_1 * denominator1_2);
			}
			if (denominator2_1 > 0 && denominator2_2 > 0) {
				basisSecondDiff += (-_degree) * (-_degree + 1) / (denominator2_1 * denominator2_2) * basisFunction(j + 1, t, _degree - 2);
			}
			if (denominator1_1 > 0 && denominator2_2 > 0) {
				basisSecondDiff += (-_degree) / (denominator1_1 * denominator2_2) * basisFunction(j, t, _degree - 1);
			}
			if (denominator2_1 > 0 && denominator1_2 > 0) {
				basisSecondDiff += _degree / (denominator2_1 * denominator1_2) * basisFunction(j + 1, t, _degree - 1);
			}
		}

		ans += basisSecondDiff * controlPoint[j];

	}

	return ans;

}

/*
*@brief: Check if the knot vector is correct.
***/
void BSplineCurve::isKnotCorrect()
{
	for (int i = 0; i < knotVector.size() - 1; i++)
	{
		assert(knotVector[i] <= knotVector[i + 1]);
	}
	Eigen::Vector2d point(0.0, 0.0);

}

/*
	B样条曲线基转换矩阵。阶数不高(节点比较多)的画可以一阶一阶的降，一直降到四阶B样条曲线
	阶数不高且要插入的节点较少可以使用节点插入算法使其成为分段Bezier曲线，再降阶
	如果阶数比较高，潘日晶——B样条曲线最小二乘降阶方法一次性降多阶
*/
