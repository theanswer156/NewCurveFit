#include "BezierCurve.h"
#include "BiArc.h"

/**
*	@brief 基于数值积分，计算曲线在指定区间上的弧长
*	@param beginT 区间起点(包含起点)
*	@param endT 区间终点(包含终点)
*   @return 曲线在指定区间上的弧长.

* 
*/
// 如果曲线包含有尖角或者变化很大的弯曲，则此方法可能计算不准确
// 但是不影响我们利用贝塞尔曲线来拟合离散的数据点
// 因为我们在拟合计算误差的时候，采用了弧长参数化方法(等距分段)，然后计算误差
// 如果有尖角或者变化很大的弯曲，最大误差点很可能是在尖角处或者变化很大的弯曲处
// 而我们采用的方法是从最大误差处进行分割，所以不会有问题
void BezierCurve::doComputeLength(const std::vector<Vector2d>& contralPoints, const double& beginT, const double& endT)
{
	const double TOLEREANCE = 1e-6;
	const double midT = (beginT + endT) / 2.0;
	this->m_length = adaptive_simpson_3_8(contralPoints,beginT, endT, TOLEREANCE);
}


/* @brief 依据建立的模型，计算长度系数
* @param contralPoints 控制点向量
* @param beginT 区间起点下标(包含起点)
* @param endT	区间终点下标(包含终点)
* @return		离散数据点的贝塞尔曲线的四个控制点
* @note			该方法尚有缺陷，需要进一步完善
*
*/
std::vector<Vector2d> BezierCurve::bezierGenerate(std::vector<Vector2d>& points_, const int& begin, const int& end,vector<double>& parameter)
{
	std::vector<Vector2d> ans(4, Vector2d::Zero());
	int size = end - begin + 1;

	if (size == 1)
	{
		return ans;
	}

	//! 记录第一个控制点和第四个控制点
	ans[0] = points_[begin];
	ans[3] = points_[end];

	//! 使用差分公式计算起始点和终止点的单位切向
	//! 这里只使用了前后想一阶差分公式当作实例以示算法的可行性
	//! 但是由于数据点的弯曲程度不同，切向的计算可能存在问题
	//! 后面只能使用更加精确的切向计算方法以保证结果的准确性
	//! 
	//! 切向的计算要不要划分成中心差分公式或者更高阶的前向或者后向差分公式
	Eigen::Vector2d leftTangent = this->tangent[begin];
	Eigen::Vector2d rightTangent = -this->tangent[end];

	if (size == 2)
	{
		double dist = (points_[end] - points_[begin]).norm()/3.0;
		ans[1] = points_[begin] + leftTangent * dist;
		ans[2] = points_[end] + rightTangent * dist;
		return ans;
	}


	double B = 0.0;
	double C = 0.0;
	double BC = 0.0;
	double AB = 0.0;
	double AC = 0.0;
	for (size_t i = 1; i < size; ++i)
	{
		double t = parameter[i];
		Eigen::Vector2d Pa = points_[i] - ans[0] - std::pow(t, 2) * (3 - 2 * t) * (ans[3] - ans[0]);
		Eigen::Vector2d Pb = -3 * t * std::pow(1 - t, 2) * leftTangent;
		Eigen::Vector2d Pc = -3 * std::pow(t, 2) * (1 - t) * rightTangent;


		B+=Pb.squaredNorm();
		C+=Pc.squaredNorm();
		BC+=Pb.dot(Pc);
		AB+=Pa.dot(Pb);
		AC+=Pa.dot(Pc);
	}
	double judge = (B * C - BC * BC);


	//! B * C - BC * BC == 0 说明点在同一条直线上，但是我们已经过滤了直线
	//! 所以这种情况不会出现在我们的计算当中
	if (judge == 0)
	{
		std::cout << "Error: B * C - BC * BC == 0." << "\t" << "begin:" << begin << "\t" << "end:" << end << std::endl;
	}

	double Ra = judge == 0 ? 0 : std::abs((AB * C - BC * AC) / judge);
	double Rb = judge == 0 ? 0 : std::abs((AC * B - AB * BC) / judge);


	//! 加入容错机制，如果Ra或者Rb大于曲线长度的2倍，则认为曲线存在异常，
	//! 则将Ra和Rb设置为曲线长度的1/3
	double dDist = (ans[3] - ans[0]).norm();
	bool bAbnormal = Ra > dDist * 2.0 || Rb > dDist * 2.0;
	double epsilon = 1e-6;
	if (bAbnormal || Ra < epsilon || Rb < epsilon)
	{
		Ra = dDist / 3.0;
		Rb = dDist / 3.0;
	}
	Eigen::Vector2d contralPoint1 = ans[0] + Ra * leftTangent;
	Eigen::Vector2d contralPoint2 = ans[3] + Rb * rightTangent;


	ans[1] = contralPoint1;
	ans[2] = contralPoint2;


	return ans;
}

/*
* @brief 计算曲线在指定区间上的弧长参数化
* @param points_	离散数据点集合
* @param begin		离散数据点起点下标(包含起点)
* @param end		离散数据点终点下标(包含终点)
* @return			弧长参数化序列
**/
vector<double> BezierCurve::arcLengthParameterize(std::vector<Vector2d>& points_, const int& begin, const int& end)
{
	assert(begin >= 0 && end < points_.size());
	//! tPrime和dist复用
	std::vector<double> tPrime{0.0};

	double dist = 0.0;
	for (int i = begin + 1; i <= end; i++)
	{
		dist+= (points_[i] - points_[i-1]).norm();
		tPrime.push_back(dist);
	}
	for (int i = 0; i < tPrime.size(); i++)
	{
		tPrime[i] /= dist;
	}
	return tPrime;
}

/*
* @brief 在第一次拟合误差稍大于容差，但是又超出不大时，
*	     尝试更新弧长参数，尽量避免使用过多的贝塞尔曲线进行拟合
* @param contralPoints	控制点向量
* @param points			参数化的离散数据点
* @param begin			离散数据点起点下标(包含起点)
* @param end			离散数据点终点下标(包含终点)
* @return parameter		引用传递重新参数化的弧长参数化序列
* @note MAX_ITERATION	最大迭代次数(没用上)
* @note MAX_STEPS		平均步长(限制最大的变化范围)
* @note					当我们再次对parameter进行参数化时，就与弦长的比没有关系了，
*						此时，我们优化的是离散数据点与曲线上的最短距离，并计算对应的时间t
**/
void BezierCurve::reParameterize(const std::vector<Eigen::Vector2d>& contralPoints, const std::vector<Eigen::Vector2d>& points,const int& begin, const int& end, std::vector<double>& parameter)
{

	const int MAX_ITERATION = 2;
	int size = end - begin + 1;
	double MAX_STEPS = 1.0 / size;
	assert(size == parameter.size());

	for (int iter = 0; iter < MAX_ITERATION; iter++)
	{
		//! 最小化每个点到贝塞尔曲线的距离，求得新的参数化序列
		for (size_t i = 0; i < parameter.size();++i)
		{
			//! 计算贝塞尔曲线上对应点和离散数据点的向量
			double tPrime = parameter[i];
			Eigen::Vector2d point = pointAt(contralPoints, tPrime);
			Eigen::Vector2d Vec = point - points[i + begin];

			Eigen::Vector2d Derivative = computeDerivative(contralPoints, tPrime);
			Eigen::Vector2d DerivativePrime = computePrimePrime(contralPoints, tPrime);


			double f = Vec.dot(Derivative);
			double fPrime = Derivative.squaredNorm()+Vec.dot(DerivativePrime);
			
			//! 计算调整后的时间T，并将变化范围限制在[-0.1,0.1]之间
			if (fPrime != 0.0)
			{
				parameter[i] -= std::clamp(f / fPrime, -MAX_STEPS, MAX_STEPS);
				parameter[i] = std::clamp(parameter[i], 0.0, 1.0);
			}
		}
	}
}

/*
* @brief 计算贝塞尔曲线与离散数据点的最大平方误差
* @param controlPoints	控制点向量
* @param points_		离散数据点集合
* @param begin			离散数据点起点下标(包含起点)
* @param end			离散数据点终点下标(包含终点)
* @param parameter		弧长参数化序列
* @return				最大误差值，并且将最大误差的相对偏移量以引用的形式返回
* @note					不能将maxError当作引用传进去，只能将maxErrorIndex使用引用传进去，
*/
double BezierCurve::computeMaxError(std::vector<Vector2d>& contralPoints, std::vector<Vector2d>& points_, const int& begin, const int& end,vector<double>& parameter,int& maxErrorIndex)
{
	double maxError = 0.0;
	for (size_t i = begin; i <= end; i++)
	{
		//! 计算贝塞尔曲线上对应点和离散数据点的向量
		Eigen::Vector2d point = pointAt(contralPoints,parameter[i-begin]);
		double error = (point - points_[i]).squaredNorm();
		if (error > maxError)
		{
			maxError = error;
			maxErrorIndex = i;
		}
	}
	
	return maxError;
}


/*
* @brief 在经过多次reparameterize误差优化以后如果误差任然很大，
*        则尝试对数据点进行平滑处理，以便得到更加精确的曲线
* @param points_		离散数据点集合
* @param begin			离散数据点起点下标(包含起点)
* @param end			离散数据点终点下标(包含终点)
* @return vector<Vector2d> 平滑后的数据
* @note					该函数后续可能跟随数据稠密一起作用，
*						因为点越少，算法效果越差。
* 
**/
vector<Vector2d> BezierCurve::smoothData(vector<Vector2d>& points_, const int& begin, const int& end)
{
	vector<Vector2d> ans;
	ans.emplace_back(points_[begin]);

	//! 使用周围三个点对数据进行平均平滑，注意：起始点和终止点不参与平滑
	for (int i = begin + 1; i < end; i++)
	{
		Eigen::Vector2d point = (points_[i - 1] + points_[i] + points_[i + 1]) / 3.0;
		ans.emplace_back(point);
	}

	ans.emplace_back(points_[end]);

	return ans;
}

BezierCurve::BezierCurve(const std::vector<Vector2d>& _points)
{	
	denseData(_points);
	doComputeTangent();
	bezierCurveFitting(points, 0, points.size() - 1);
}

BezierCurve::BezierCurve(const std::vector<Vector2d>& _points,const std::vector<Vector2d>& _tangent)
{
	this->points = _points;
	this->tangent = _tangent;
	//denseData(_points);
	//doComputeTangent();
	bezierCurveFitting(points, 0, points.size() - 1);
}

/*
* @brief 使用贝塞尔曲线拟合离散数据点
* @param points_		离散数据点集合
* @param begin			离散数据点起点下标(包含起点)
* @param end			离散数据点终点下标(包含终点)
* @param TOLERENCE		容差
* @param MAX_ITERATION	最大迭代次数
* @return vector<Vector2d> 拟合整条曲线的所有的贝塞尔曲线组成的控制点向量
* @note 该方法尚未完成，需要进一步完善
* @note 点太少了又不能发挥正确的作用  点太少了就用直线？？？
***/
void BezierCurve::bezierCurveFitting(std::vector<Vector2d>& points_, const int& begin, const int& end, const double& TOLERENCE)
{
	//! 设置返回向量
	std::vector<Eigen::Vector2d> ans(4);
	const int MAX_ITERATION = 10;
	int size = end - begin + 1;
	if (size == 2)
	{

		ans[0]=(points_[begin]);
		ans[1]=(points_[begin]);
		ans[2]=(points_[end]);
		ans[3]=(points_[end]);
		this->ctrlPoints.insert(this->ctrlPoints.end(), ans.begin(), ans.end());
		return;
	}

	//! 计算弧长参数化序列,第一次参数化、第一次拟合，所以只需要离散数据点集
	std::vector<double> tPrime = arcLengthParameterize(points_, begin, end);
	std::vector<Eigen::Vector2d> controlPoints = bezierGenerate(points_, begin, end,tPrime);

	//! 计算最大误差，并使用引用记录最大误差点的相对偏移量
	int maxErrorIndex = -1;
	double maxError = computeMaxError(controlPoints, points_, begin, end, tPrime,maxErrorIndex);


	//! 如果误差在容差范围内，则直接返回控制点向量
	//! 从而使用一条贝塞尔曲线就能够完全拟合离散数据点
	if (maxError < TOLERENCE)
	{
		this->ctrlPoints.insert(this->ctrlPoints.end(), controlPoints.begin(), controlPoints.end());
		return;
	}

	//! 如果误差很大，超过设置容差的一定范围，那么直接从最大误差处分割数据点
	//! 然后再次进行拟合，直到误差小于容差为止
	if (maxError > 50.0 * TOLERENCE)
	{
		bezierCurveFitting(points_, begin, maxErrorIndex, TOLERENCE);
		bezierCurveFitting(points_, maxErrorIndex, end, TOLERENCE);
		return;
	}


	//! 如果误差在容差的10倍范围内，由于容差设置的也很小，
	//! 所以我们希望能够直接进行优化以便使用尽量少的贝塞尔曲线拟合数据点
	//! 如果调整后能够到能接受的误差范围,则直接返回控制点向量  这里我们设置为2倍容差
	//if (maxError <= 10.0 * TOLERENCE)
	//{
	//}
	for (int i = 0; i < MAX_ITERATION; ++i)
	{
		reParameterize(controlPoints, points_, begin, end, tPrime);
		controlPoints = bezierGenerate(points_, begin, end,tPrime);
		double max_error = computeMaxError(controlPoints, points_, begin, end, tPrime,maxErrorIndex);
		if (max_error < 2.0 * TOLERENCE)
		{
			this->ctrlPoints.insert(this->ctrlPoints.end(), controlPoints.begin(), controlPoints.end());
			return;
		}
	}


	bezierCurveFitting(points_, begin, maxErrorIndex, TOLERENCE);
	bezierCurveFitting(points_, maxErrorIndex, end, TOLERENCE);

//! 如果经过误差调整后还是很大，试着对数据点进行平滑处理，然后再次进行拟合
//! 若是平滑后拟合仍然不理想，则将数据点分割成两部分，分别进行拟合
//! 
//! 从结果来看，平滑后可以稍微解决问题，但是仍然存在一些问题。
//! 例如，我们可能会失去曲线的一些特性，例如弯曲程度、尖角等。


//vector<Vector2d> smoothedPoints = smoothData(points_, begin, end);
//
//if (!bezierCurveFittingWithSmooth(smoothedPoints, 0, smoothedPoints.size() - 1, TOLERENCE))
//{
//	//! 如果误差超出容差的10倍范围，我们通过之前计算的最大误差点的偏移量
//	//! 将数据点分割成两部分，分别进行拟合
//	bezierCurveFitting(points_, begin, maxErrorIndex, TOLERENCE);
//	bezierCurveFitting(points_, maxErrorIndex, end, TOLERENCE);
//}

}

bool BezierCurve::bezierCurveFittingWithSmooth(std::vector<Vector2d>& points_, const int& begin, const int& end, const double& TOLERENCE)
{
	//! 设置返回向量
	std::vector<Eigen::Vector2d> ans(4);
	const int MAX_ITERATION = 10;
	int size = end - begin + 1;
	if (size == 2)
	{

		ans.emplace_back(points_[begin]);
		ans.emplace_back(points_[begin + 1]);
		ans.emplace_back(points_[end]);
		ans.emplace_back(points_[end - 1]);
		this->ctrlPoints.insert(this->ctrlPoints.end(), ans.begin(), ans.end());
		return true;
	}

	//! 计算弧长参数化序列,第一次参数化、第一次拟合，所以只需要离散数据点集
	std::vector<double> tPrime = arcLengthParameterize(points_, begin, end);
	std::vector<Eigen::Vector2d> controlPoints = bezierGenerate(points_, begin, end, tPrime);

	//! 计算最大误差，并使用引用记录最大误差点的相对偏移量
	int maxErrorIndex = -1;
	double maxError = computeMaxError(controlPoints, points_, begin, end, tPrime, maxErrorIndex);


	//! 如果误差在容差范围内，则直接返回控制点向量
	//! 从而使用一条贝塞尔曲线就能够完全拟合离散数据点
	if (maxError < TOLERENCE)
	{
		this->ctrlPoints.insert(this->ctrlPoints.end(), controlPoints.begin(), controlPoints.end());
		return true;
	}


	//! 如果误差在容差的10倍范围内，由于容差设置的也很小，
	//! 所以我们希望能够直接进行优化以便使用尽量少的贝塞尔曲线拟合数据点
	//! 如果调整后能够到能接受的误差范围,则直接返回控制点向量  这里我们设置为2倍容差
	//if (maxError <= 10.0 * TOLERENCE)
	//{
	//}
	for (int i = 0; i < MAX_ITERATION; ++i)
	{
		reParameterize(controlPoints, points_, begin, end, tPrime);
		controlPoints = bezierGenerate(points_, begin, end, tPrime);
		double max_error = computeMaxError(controlPoints, points_, begin, end, tPrime, maxErrorIndex);
		if (max_error < 2.0 * TOLERENCE)
		{
			this->ctrlPoints.insert(this->ctrlPoints.end(), controlPoints.begin(), controlPoints.end());
			return true;
		}
	}


	return false;
}

/*
*  @brief 将符合计算结果的控制点输出到ctrPoints中
*  @param ctrPoints 控制点向量
*/
std::vector<Vector2d> BezierCurve::outContralPoints() const
{
	return this->ctrlPoints;
}

/*
*@brief 设置控制点
*@param ctrlPoints 控制点向量
***/
void BezierCurve::setControlPoints(const std::vector<Vector2d>& ctrlPoints)
{
	assert(ctrlPoints.size() == 4);
	this->ctrlPoints = ctrlPoints;
}

/*
* @brief 向距离较远的两个点之中添加新的点
****/
void BezierCurve::denseData(const vector<Vector2d>& points_)
{
	assert(!points_.empty());
	size_t n = points_.size();
	int iBegin = 0;
	int iEnd = 1;
	const double MIN_DISTANCE = 2.0;
	while (iEnd < n)
	{
		double dist = (points_[iEnd] - points_[iBegin]).norm();
		if (dist != 0.0)
		{
			int INSERTNUMBER = dist < MIN_DISTANCE ? 0 : static_cast<int>(std::ceil(dist / MIN_DISTANCE));;
			this->points.emplace_back(points_[iBegin]);
			Eigen::Vector2d direc = (points_[iEnd] - points_[iBegin]);
			for (int i = 1; i < INSERTNUMBER; i++)
			{
				Eigen::Vector2d insertPoint = points_[iBegin] + (i * 1.0 / INSERTNUMBER) * direc;
				this->points.emplace_back(insertPoint);
			}
		}
		++iEnd;
		++iBegin;
	}
	std::cout << string(50,'-') << std::endl;
	std::cout << "denseData size:" << this->points.size() << std::endl;
	std::cout << string(50, '-') << std::endl;
}

/*
* @brief 等间距数值方法计算所有点的切向量
*
***/
void BezierCurve::doComputeTangent()
{
	assert(!this->points.empty());
	size_t n = this->points.size();
	tangent.clear();
	Eigen::Vector2d tangent(0.0, 0.0);

	for (size_t i = 0; i < n; i++)
	{
		if (i > 0 && i < n - 1)
		{
			tangent = (this->points[i + 1] - this->points[i - 1]) / 2.0;
		}
		else if (i == 0)
		{
			tangent = this->points[1] - this->points[0];
		}
		else
		{
			tangent = this->points[n-1] - this->points[n-2];
		}
		tangent.normalize();
		this->tangent.emplace_back(tangent);
	}
}

/*
* @brief 使用增量式最小二乘法计算切向
****/
void BezierCurve::doComputeTangent_LSE()
{
	assert(!this->points.empty() && !this->tangent.empty());
	double a = 0;
	double b = 0;
	double X = 0;
	double Y = 0;
	double XY = 0;
	double XX = 0;
	int MINLENGTH = 7;
	int WINDOWSIZE = 9;
	int HALFWINDOWSIZE = std::floor(WINDOWSIZE / 2);
	size_t iBegin = 0;
	size_t iEnd = static_cast<size_t> (WINDOWSIZE - 1);
	size_t iMid = static_cast<size_t> (std::floor(iBegin + iEnd) / 2);

	for (size_t index = 0; index < iMid; index++)
	{
		this->tangent[index] = (this->points[index+1] - this->points[index]).normalized();
	}

	for (size_t index = 0; index < WINDOWSIZE; index++)
	{
		X += this->points[index].x();
		Y += this->points[index].y();
		XY += this->points[index].x() * this->points[index].y();
		XX += this->points[index].x() * this->points[index].x();
	}

	while (iEnd < this->points.size())
	{
		bool bAbNormal = false;
		if (iBegin)
		{
			X -= (this->points[iBegin - 1].x());
			X += (this->points[iEnd].x());
			Y -= (this->points[iBegin - 1].y());
			Y += (this->points[iEnd].y());
			XX -= (this->points[iBegin - 1].x() * this->points[iBegin - 1].x());
			XX += (this->points[iEnd].x() * this->points[iEnd].x());
			XY -= (this->points[iBegin - 1].x() * this->points[iBegin - 1].y());
			XY += (this->points[iEnd].x() * this->points[iEnd].y());
		}


		if (std::abs(WINDOWSIZE * XX - X * X) < 1e-4 || std::abs(WINDOWSIZE * XY - X * Y) < 1e-4)
		{
			std::cout << "IMID : " << iMid << "\t" << "IBEGIN : " << iBegin << "\t" << "IEND : " << iEnd << std::endl;
			bAbNormal = true;
			a = 0.0;
			b = 0.0;
		}
		else
		{
			a = (WINDOWSIZE * XY - X * Y) / (WINDOWSIZE * XX - X * X);
			b = (Y - a * X) / (WINDOWSIZE);
		}


		//!	根据拟合的直线的斜率得到该点的单位切向量
		//! 注意切线的方向与曲线的方向相同
		{
			Eigen::Vector2d lineVector = this->points[iEnd] - this->points[iMid];
			if (bAbNormal)
			{
				bool bXDirection = std::abs(lineVector.x()) > std::abs(lineVector.y());
				double dValue = 0.0;
				if (bXDirection)
				{
					dValue = lineVector.x() > 0 ? 1 : -1;
					this->tangent[iMid] = Eigen::Vector2d{ dValue,0 };
				}
				else
				{
					dValue = lineVector.y() > 0 ? 1 : -1;
					this->tangent[iMid] = Eigen::Vector2d{ 0, dValue };
				}
			}
			else
			{
				//!	默认Y的分量为正
				double modulus = std::sqrt(1 + a * a);
				bool bSameDirect = lineVector.x() * a + lineVector.y() >= 0;
				//this->tangent[iMid] = bSameDirect ? Point{ 1 / modulus,a / modulus } : Point{ -1 / modulus,-a / modulus };
				this->tangent[iMid] = Eigen::Vector2d{ 1 / modulus,a / modulus };
				double angle = std::atan2(this->tangent[iMid].x() * lineVector.y() - this->tangent[iMid].y() * lineVector.x(), this->tangent[iMid].dot(lineVector));
				this->tangent[iMid] = std::abs(angle) <= PI / 2 ? this->tangent[iMid] : -this->tangent[iMid];

			}
		}
	}
}

/*
*  @brief 计算曲线在指定点处的1阶导数
*  @param contralPoints 控制点向量
*  @param t ([0,1]之间)指定点的时间
*  @return 曲线在指定点处的1阶导数
*/
inline Eigen::Vector2d BezierCurve::computeDerivative(const std::vector<Eigen::Vector2d>& contralPoints, const double& t)
{
	Eigen::Vector2d result;
	assert(!contralPoints.empty());
	double newT = std::clamp(t, 0.0, 1.0);
	result = 3 * std::pow(1 - newT, 2) * (contralPoints[1] - contralPoints[0])
			+ 6 * (1 - newT) * newT * (contralPoints[2] - contralPoints[1])
			+ 3 * std::pow(newT, 2) * (contralPoints[3] - contralPoints[2]);

	return result;
}

/*
*  @brief 计算曲线在指定点处的2阶导数
*  @param contralPoints 控制点向量
*  @param t ([0,1]之间)指定点的时间
*  @return 曲线在指定点处的2阶导数
*/
inline Eigen::Vector2d BezierCurve::computePrimePrime(const std::vector<Eigen::Vector2d>& contralPoints, const double& t)
{
	Eigen::Vector2d result;
	double newT = std::clamp(t, 0.0, 1.0);
	result = 6 * (1 - newT) * (contralPoints[2] - 2 * contralPoints[1] + contralPoints[0])
		   + 6 * newT * (contralPoints[3] - 2 * contralPoints[2] + contralPoints[1]);

	return result;
}
/*
* @brief 计算曲线在指定点处的位置
* @param controlPoints 控制点向量
* @param t ([0,1]之间)指定点的时间
* @return 曲线在指定点处的位置
*/
inline Eigen::Vector2d BezierCurve::pointAt(const std::vector<Vector2d>& controlPoints,const double& t)
{
	assert(t >= 0.0 && t <= 1.0);
	if (t == 0.0)
	{
		return controlPoints[0];
	}
	if (t == 1.0)
	{
		return controlPoints[3];
	}
	Eigen::Vector2d result;

	//double newT = std::clamp(t, 0.0, 1.0);

	result = controlPoints[0] * std::pow(1 - t, 3) + controlPoints[1] * 3 * std::pow(1 - t, 2) * t 
		   + controlPoints[2] * 3 * (1 - t) * std::pow(t, 2) + controlPoints[3] * std::pow(t, 3);
	return result;
}


/*
*  @brief 采用辛普森3/8方法计算曲线在指定区间上的弧长
*  @param leftTime 区间起点(包含起点)
*  @param rightTime 区间终点(包含终点)
*  @return 曲线在指定区间上的弧长
*/
double BezierCurve::simpson_3_8(const std::vector<Vector2d>& contralPoints,const double& leftTime, const double& rightTime)
{
	double midTime_L = (2 * leftTime + rightTime) / 3.0;
	double midTime_R = (leftTime + 2 * rightTime) / 3.0;
	auto p1 = computeDerivative(contralPoints,leftTime).norm();
	auto p2 = 3 * computeDerivative(contralPoints, midTime_L).norm();
	auto p3 = 3 * computeDerivative(contralPoints, midTime_R).norm();
	auto p4 = computeDerivative(contralPoints, rightTime).norm();


	return (p1+p2+p3+p4)*(rightTime-leftTime)/8.0;
}

double BezierCurve::adaptive_simpson_3_8(const std::vector<Vector2d>& contralPoints, const double& leftTime, const double& rightTime, const double& TOLERENCE)
{
	const double midTime = (leftTime + rightTime) / 2.0;
	double simpson_Total = simpson_3_8(contralPoints, leftTime, rightTime);
	double simpson_Left = simpson_3_8(contralPoints, leftTime, midTime);
	double simpson_Right = simpson_3_8(contralPoints, midTime, rightTime);
	double error = simpson_Left + simpson_Right - simpson_Total;
	if (std::abs(error) <= 15.0 * TOLERENCE)
	{
		return simpson_Left + simpson_Right + error / 15.0;
	}
	return adaptive_simpson_3_8(contralPoints,leftTime, midTime, TOLERENCE / 2.0)
		  + adaptive_simpson_3_8(contralPoints,midTime, rightTime, TOLERENCE / 2.0);
}
