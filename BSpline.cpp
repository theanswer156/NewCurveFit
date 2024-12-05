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


/*
*@berif: Boehm算法插入节点，每次只能插入一个节点，计算控制点
*		 然后把插入的节点和插入到初始节点向量，计算出来的控制点放到对应位置
*		 然后再一次使用Boehm算法进行下一个节点的插入
***/
void BSplineCurve::BoehmKnotInsert()
{
	//! 确定B样条的有效域
	int beginKnotIndex = this->degree - 1;
	int endKnotIndex = this->controlPoint.size();
	
	//! 确定插入节点索引以及重复度
	unordered_map<double, int> knotIndexMap;
	//! 先试一下所有的节点都插入  这很符合贝塞尔曲线的节点的重复度
	//! 这里还真不能统计所有的节点，只能在B样条有效域内统计
	for (int i = beginKnotIndex; i <= endKnotIndex; i++)
	{
		knotIndexMap[knotVector[i]]++;
	}

	vector<pair<double, int>> kontInsert(knotIndexMap.begin(), knotIndexMap.end());


	for (const auto& _pair : kontInsert)
	{
		double knot = _pair.first;
		int repeat = _pair.second;
		//! 插入次数和插入位置
		int insertTimes = this->degree - repeat;
		//! 这里插入位置是找到第一个大于等于knot的位置
		auto insertPos = std::distance(lower_bound(this->knotVector.begin(), this->knotVector.end(), knot),this->knotVector.begin());

		assert(insertTimes >= 0 && "Insert times must be greater than or equal to 0");

		vector<Vector2d> subControlPoints(this->controlPoint.begin() + insertPos - this->degree, controlPoint.begin() + insertPos + 2);
		vector<double> subKnotVector(this->knotVector.begin() + insertPos - this->degree - 1, this->knotVector.begin() + insertPos + this->degree + 2);

		auto subControlPoints_size = subControlPoints.size();
		auto subKnotVector_size = subKnotVector.size();

		//! 一次性插入insertTimes次节点knot  我只要申请一个subInsertControlPoints的空间
		//! 以此存储插入后的控制点就行了  对于K次多项式，如果插入节点的重复度为r，则需要重新计算K+r-1个控制点
		//! 因此这里只需要申请insertTimes+degree个控制点的空间
		//! 此时subControlPoints的长度为degree+1，我们可以在subControlPoints的基础上进行开花
		//! 每次都正好从layer开始，将计算出来的控制点就此放入，如果往后面放入就覆盖了后面要计算的控制点
		vector<Vector2d> subInsertControlPoints(insertTimes+this->degree, Vector2d(0.0, 0.0));
		int subInsertControlPoints_size = subInsertControlPoints.size();

		for (int layer = 1; layer <= insertTimes; layer++)
		{
			//! 根据开花算法以此计算每一层的控制点
			for (int i = 0; i <= subControlPoints_size - layer; i++)
			{
				double normalizeCoeff = (subKnotVector[this->degree + i] - subKnotVector[layer - 1 + i]);
				subControlPoints[i] = (subKnotVector[this->degree + i] - knot) * subControlPoints[i];
				subControlPoints[i]+= (knot - subKnotVector[layer - 1 + i])*subControlPoints[i + 1];
				subControlPoints[i] /= normalizeCoeff;
			}

			//! 如果节点插入到了最后一层，将最后一层所有的控制点都放入subInsertControlPoints中
			//! 否则，将第一个计算和最后一个更新的点放入subInsertControlPoints中
			if (layer == insertTimes)
			{
				int endIndex = subControlPoints_size - layer - 1;
				for (int j = 0; j <= endIndex; j++)
				{
					subInsertControlPoints[j + layer - 1] = subControlPoints[j];
				}
			}
			else
			{
				subInsertControlPoints[layer-1] = subControlPoints[0];
				subInsertControlPoints[subInsertControlPoints_size - layer - 1] = subControlPoints[subControlPoints_size - layer - 1];

			}
		}
		//! 最后将插入后的控制点放入原来的控制点中  但是频繁的往中间插入效率不高，还是修改成为每次往最后插入节点
		
		
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

/*
*@brief: B样条曲线基函数
* @param i: 控制点序号
* @param t: 时间参数
* @param _degree: 基函数次数(不是阶数)
***/
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
