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
	double startT = *(knotVector.begin() + degree - 1);
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
	//! 这里还真不能统计所有的节点，只能在B样条有效域内统计 特别的
	//! 这个算法仅仅适用于准均匀B样条曲线  下面要做的就是能够找一个
	//! 可以在有效域外进行节点插入的算法或者是找一个算法使得
	//! 能够将任意B样条曲线变成准均匀的B样条曲线
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
		int insertTimes = this->degree - repeat - 1;
		//! 这里插入位置是找到第一个大于等于knot的位置
		auto insertPos = std::distance(lower_bound(this->knotVector.begin(), this->knotVector.end(), knot), this->knotVector.begin());

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
		vector<Vector2d> subInsertControlPoints(insertTimes + this->degree, Vector2d(0.0, 0.0));
		int subInsertControlPoints_size = subInsertControlPoints.size();

		for (int layer = 1; layer <= insertTimes; layer++)
		{
			//! 根据开花算法以此计算每一层的控制点
			for (int i = 0; i <= subControlPoints_size - layer; i++)
			{
				double normalizeCoeff = (subKnotVector[this->degree + i] - subKnotVector[layer - 1 + i]);
				if (normalizeCoeff != 0.0)
				{
					subControlPoints[i] = (subKnotVector[this->degree + i] - knot) * subControlPoints[i];
					subControlPoints[i] += (knot - subKnotVector[layer - 1 + i]) * subControlPoints[i + 1];
					subControlPoints[i] /= normalizeCoeff;
				}
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
				subInsertControlPoints[layer - 1] = subControlPoints[0];
				subInsertControlPoints[subInsertControlPoints_size - layer - 1] = subControlPoints[subControlPoints_size - layer - 1];

			}
		}
		//! 最后将插入后的控制点放入原来的控制点中  但是频繁的往中间插入效率不高，还是修改成为每次往最后插入节点


	}



}

/*
*@brief: Qin算法插入节点，每次只能插入一个节点，计算所有控制点
*		 将输入的样条曲线通过Qin节点插入算法转化为准均匀的样条曲线
*		 然后使用boehm节点插入算法转化为分段Bezier曲线
*
***/
BSplineCurve BSplineCurve::QinKnotInsert(const BSplineCurve& curve, const double& insertKnot)
{
	//! 首先确定节点插入前后B样条基函数的转换矩阵
	//! 
	size_t n = curve.controlPoint.size() - 1;
	int k = curve.knotVector.size() - curve.controlPoint.size() - 1;
	std::cout << "controlPoint size: " << curve.controlPoint.size() << std::endl;
	std::cout << "knotVector size: " << curve.knotVector.size() << std::endl;
	std::cout << "k: " << k << std::endl;
	Eigen::MatrixXd transMat = Eigen::MatrixXd::Zero(n + 1, n + 1);
	Eigen::MatrixXd pointMat(n + 1, 2);
	for (size_t i = 0; i < pointMat.rows(); i++)
	{
		pointMat.row(i) = curve.controlPoint[i].transpose();
	}
	for (size_t j = 0; j < transMat.rows(); j++)
	{
		double Alpha = 0.0;
		double Beta = 0.0;
		double delta_1 = (curve.knotVector[j + k - 1] - curve.knotVector[j]);
		double delta_2 = (curve.knotVector[j + k] - curve.knotVector[j + 1]);
		if (insertKnot <= curve.knotVector[j])
		{
			Beta = 1.0;
		}
		else if (insertKnot <= curve.knotVector[j + 1])
		{
			Alpha = delta_1 == 0 ? 0 : (insertKnot - curve.knotVector[j]) / delta_1;
			Beta = 1.0;
		}
		else if (insertKnot <= curve.knotVector[j + k - 1])
		{
			Alpha = delta_1 == 0 ? 0 : (insertKnot - curve.knotVector[j]) / delta_1;
			Beta = delta_2 == 0 ? 0 : (curve.knotVector[j + k] - insertKnot) / delta_2;
		}
		else if (insertKnot <= curve.knotVector[j + k])
		{
			Alpha = 1.0;
			Beta = delta_2 == 0 ? 0 : (curve.knotVector[j + k] - insertKnot) / delta_2;
		}
		else
		{
			Alpha = 1.0;
		}
		transMat(j, j) = Alpha;
		if (j < n)
		{
			transMat(j, j + 1) = Beta;
		}
	}

	std::cout << "transMat:\n" << transMat << std::endl;

	BSplineCurve ans;
	//auto L = std::distance(curve.knotVector.begin(), curve.knotVector.end());
	auto L = std::distance(curve.knotVector.begin(), std::upper_bound(curve.knotVector.begin(), curve.knotVector.end(), insertKnot));
	if (L > (curve.knotVector.size() >> 1))
	{
		L = std::distance(curve.knotVector.begin(), std::lower_bound(curve.knotVector.begin(), curve.knotVector.end(), insertKnot));
	}
	//L > (curve.knotVector.size() >> 1) ? L -= 1 : L += 1;
	std::cout << "L: " << L << std::endl;
	if (L == 0)
	{
		return curve;
	}
	Eigen::MatrixXd subTransMat = Eigen::MatrixXd::Identity(n + 1, n + 1);
	bool isKnotIncreasing = false;

	if (insertKnot == curve.knotVector[k - 1] || L < k)
	{
		Eigen::MatrixXd mat = transMat.block(0, 1, k, k - 1).transpose();
		std::cout << "mat:\n" << mat << std::endl;
		subTransMat.block(0, 0, k - 1, k) = transMat.block(0, 1, k, k - 1).transpose();
		ans.knotVector.assign(curve.knotVector.begin() + 1, curve.knotVector.end());
	}
	else if (insertKnot == curve.knotVector[n + 1] || L > n)
	{
		Eigen::MatrixXd mat = transMat.block(n - k + 1, n - k + 2, k, k - 1).transpose();
		std::cout << "mat:\n" << mat << std::endl;
		subTransMat.block(n + 2 - k, n + 2 - k, k - 2, k - 1) = transMat.block(n - k + 1, n - k + 2, k - 1, k - 2).transpose();
		ans.knotVector.assign(curve.knotVector.begin(), curve.knotVector.end() - 1);
	}
	else
	{
		ans.knotVector.assign(curve.knotVector.begin(), curve.knotVector.end());
		isKnotIncreasing = true;
	}

	std::cout << "subTransMat:\n" << subTransMat << std::endl;

	//! 如果不属于情况①和③，节点向量增加一个，相应的控制点也增加一个
	if (isKnotIncreasing)
	{
		for (size_t j = 0; j <= L - k + 1; ++j)
		{
			ans.controlPoint.emplace_back(curve.controlPoint[j]);
		}
		for (size_t j = L - k + 2; j <= L; j++)
		{
			double alpha_j = (insertKnot - curve.knotVector[j]) / (curve.knotVector[j + k - 1] - curve.knotVector[j]);
			Eigen::Vector2d newControlPoint = (1 - alpha_j) * curve.controlPoint[j - 1] + alpha_j * curve.controlPoint[j];
			ans.controlPoint.emplace_back(newControlPoint);
		}
		for (size_t j = L + 1; j <= n; j++)
		{
			ans.controlPoint.emplace_back(curve.controlPoint[j]);
		}
	}
	else
	{
		Eigen::MatrixXd mat = subTransMat * pointMat;
		for (size_t i = 0; i < mat.rows(); i++)
		{
			ans.controlPoint.emplace_back(mat.row(i));
		}
	}
	ans.knotVector.emplace_back(insertKnot);
	std::sort(ans.knotVector.begin(), ans.knotVector.end());


	return ans;
}

void BSplineCurve::OlsoKnotInsert()
{
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
