#include "BezierToBiArc.h"
#include "BiArc.h"
#include <iostream>

/*
* @brief 构造函数
* @param controlPoints	贝塞尔曲线的控制点,四个为一组，表示一条三阶贝塞尔曲线
***/
BezierToBiArc::BezierToBiArc(const std::vector<Vector2d>& controlPoints)
{
	int n = controlPoints.size();
	assert(n % 4 == 0);
	for (size_t i = 0; i < n; i += 4)
	{
		m_vecControlPoints.emplace_back(vector<Vector2d>{controlPoints[i], controlPoints[i + 1], 
													controlPoints[i + 2], controlPoints[i + 3]});
	}
}

std::vector<BiArc> BezierToBiArc::outBiArcs() const
{
	return m_vecBiArc;
	// TODO: 在此处插入 return 语句
}

/*
* @brief 计算指定时间点处的贝塞尔曲线上的点
* @param controlPoints	贝塞尔曲线的控制点
* @param t				指定时间点
* @return				指定时间点处的贝塞尔曲线上的点
***/
inline Vector2d BezierToBiArc::PointAt(const std::vector<Vector2d>& controlPoints, const double& t)
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

	result = std::pow(1 - t, 3) * controlPoints[0] + 3 * std::pow(1 - t, 2) * t * controlPoints[1]
		+ 3 * (1 - t) * std::pow(t, 2) * controlPoints[2] + std::pow(t, 3) * controlPoints[3];
	return result;
}


/*
* @brief 计算指定时间点处的贝塞尔曲线的切线
* @param controlPoints	贝塞尔曲线的控制点
* @param t				指定时间点
* @return 指定时间点处的贝塞尔曲线的切线
***/
inline Eigen::Vector2d BezierToBiArc::computeDerivative(const vector<Vector2d>& contralPoints, const double& t)
{
	Eigen::Vector2d result;
	result = 3 * std::pow(1 - t, 2) * (contralPoints[1] - contralPoints[0])
		+ 6 * (1 - t) * t * (contralPoints[2] - contralPoints[1])
		+ 3 * std::pow(t, 2) * (contralPoints[3] - contralPoints[2]);

	return result;
}


/*
* @brief 计算三角形的内心
* @param p1	三角形的第一个顶点
* @param p2	三角形的第二个顶点
* @param p3	三角形的第三个顶点
* @return		三角形的内心
* @note		注意要判断三角形是否共线，如果共线，则返回一个空向量或者异常
*           因为几乎共线的情况是病态的，对数据异常敏感
****/
Eigen::Vector2d BezierToBiArc::innerPoint(const Vector2d& p1, const Vector2d& p2, const Vector2d& p3)
{
	Vector2d result;
	static const double INNERPOINT_TOLERENCE = 1.0-1e-5;
	//! 首先判断三角形是否共线
	Vector2d P1P2 = (p2 - p1).normalized();
	Vector2d P1P3 = (p3 - p1).normalized();
	Vector2d P2P3 = (p3 - p2).normalized();
	double d1 = P1P2.dot(P1P3);
	double d2 = P1P2.dot(P2P3);
	double d3 = P1P3.dot(P2P3);

	if (std::abs(d1) > INNERPOINT_TOLERENCE || std::abs(d2) > INNERPOINT_TOLERENCE || std::abs(d3) > INNERPOINT_TOLERENCE)
	{
		std::cout << "Error: The triangle is collinear!,we can not compute the inner point!" << std::endl;
		return Eigen::Vector2d();
	}

	//! 利用三角形内心的几何性质来计算内心坐标
	double dA = (p3-p2).norm();
	double dB = (p1-p3).norm();
	double dC = (p2-p1).norm();

	result = (dA * p1 + dB * p2 + dC * p3) / (dA + dB + dC);

	return result;
}


/*
* @brief 计算两条直线的交点
* @param p1	直线1的起点
* @param t1	直线1的方向向量(默认是单位向量)
* @param p2	直线2的起点
* @param t2	直线2的方向向量(默认是单位向量)
* @return	直线1和直线2的交点，如果两直线不相交，那么返回一个空向量或者异常
* @note		注意判断两条直线斜率是否存在、直线是否相交
*           
***/
Eigen::Vector2d BezierToBiArc::lineIntersection(const Vector2d& p1, const Vector2d& t1, const Vector2d& p2, const Vector2d& t2)
{
	static const double LINE_INTERSECTION_TOLERENCE = 1-1e-5;
	Eigen::Vector2d result;
	//! 利用两条直线的方向向量的点乘是否接近于1来判断两条直线是否相交
	//! 如果异常可以给一点perturbation，或者进行拟合的时候进行二分
	if (std::abs(t1.dot(t2)) > LINE_INTERSECTION_TOLERENCE)
	{
		std::cout << "Error: The two lines are parallel or skew!" << std::endl;
		return Eigen::Vector2d();
	}

	//! 没有异常就说明两条直线有交点，下面分别对两条直线斜率不存在的情况进行处理
	if (std::abs(t1.x()) < 1e-5)
	{
		double k = t2.y() / t2.x();
		double y = k*(p1.x() - p2.x())+p2.y();
		result = Vector2d(p1.x(), y);
	}
	else if (std::abs(t2.x()) < 1e-5)
	{
		double k = t1.y() / t1.x();
		double y = k*(p2.x() - p1.x())+p1.y();
		result = Vector2d(p2.x(), y);
	}
	else
	{
		double k1 = t1.y() / t1.x();
		double k2 = t2.y() / t2.x();
		double x = (p2.y() - p1.y() + k1 * p1.x() - k2 * p2.x()) / (k1 - k2);
		double y = k1 * (x - p1.x()) + p1.y();
		result = Vector2d(x, y);
	}

	return result;
}


/*
* @brief 计算圆弧的圆心
* @param p1	圆弧的起点
* @param t	圆弧起始点的方向向量
* @param p2	圆弧的终点
* @return	圆弧的圆心
****/
Eigen::Vector2d BezierToBiArc::computeCircleCenter(const Vector2d& p1, const Vector2d& t, const Vector2d& p2)
{
	return Eigen::Vector2d();
}

/*
* @brief	计算贝塞尔曲线曲率突变点，并以此将其分割为多个贝塞尔曲线
* @return	将分割后的贝塞尔曲线控制点存放到m_vecBezierControlPoints中
* @note		如果贝塞尔曲线有两个曲率突变点 t1,t2(t1<t2),我们从t1开始分割，
			将t1处的贝塞尔曲线分割为两个贝塞尔曲线，将第二条贝塞尔曲线再次进行分割就好了

			但是分割后的贝塞尔曲线由于局部的变化，最终可能会另外生成其他的曲率突变点。
***/
void BezierToBiArc::splitBezierByCurvature(const std::vector<Vector2d>& controlPoints)
{
	assert(controlPoints.size() == 4);
	double t1 = getCurvateMutations(controlPoints);


	if (t1 != 2.0)
	{
		Vector2d newContralPoint_0 = controlPoints[0];
		Vector2d newContralPoint_1 = (1 - t1) * controlPoints[0] + t1 * controlPoints[1];
		Vector2d newContralPoint_2 = std::pow(1 - t1, 2) * controlPoints[0] + 2 * (1 - t1) * t1 * controlPoints[1];
		Vector2d newContralPoint_3 = PointAt(controlPoints, t1);
		Vector2d newContralPoint_4 = (1 - t1) * controlPoints[2] + t1 * controlPoints[3];
		Vector2d newContralPoint_5 = std::pow(1 - t1, 2) * controlPoints[1] + 2 * (1 - t1) * t1 * controlPoints[2];
		Vector2d newContralPoint_6 = controlPoints[3];

		std::vector<Vector2d> leftControlPoints{ newContralPoint_0, newContralPoint_1, newContralPoint_2, newContralPoint_3 };
		std::vector<Vector2d> rightControlPoints{ newContralPoint_3,newContralPoint_4, newContralPoint_5, newContralPoint_6 };
		
		m_vecControlPoints.emplace_back(leftControlPoints);
		splitBezierByCurvature(rightControlPoints);
	}
	else
	{
		m_vecControlPoints.emplace_back(controlPoints);
	}

}


/*
* @brief 逐个用双圆弧拟合类中贝塞尔曲线向量，将结果存入m_vecBiArc中
* @note		这样的化方便我们在这里使用多线程加快计算速度
***/
void BezierToBiArc::fromCurvesToBiArc()
{
	for (const auto& controlPoints : m_vecControlPoints)
	{
		fromBezierToBiArc_IP(controlPoints);
	}
}

/*
* @brief 使用双圆弧拟合指定时间段内的贝塞尔曲线
* @param controlPoints	贝塞尔曲线的控制点
* @param beginTime		起点时间(默认为0.0)
* @param endTime		终点时间(默认为1.0)
* @param TOLERENCE		容差，默认为1.0
* @param POINTS_NUM		在拟合的双圆弧上POINTS_NUM个点计算拟合的误差，默认为100
* @param MIN_LENGTH		拟合圆弧的最小长度，默认为1e2，超过这个长度的双圆弧用更多的点计算误差
* @note					这里我们不方便使用之前一样的parameterize方法
***/
void BezierToBiArc::singleBezierToBiArc(const std::vector<Vector2d>& controlPoints, const double& beginTime, const double& endTime,const double& TOLERENCE)
{
	assert(controlPoints.size() == 4 && endTime > beginTime);

	double sweepTime = endTime - beginTime;
	const int POINTS_NUM = 101;
	const double MIN_LENGTH = 1e2;

	vector<double> timeSeq(POINTS_NUM, 0);
	for (int i = 0; i < POINTS_NUM; i++)
	{
		timeSeq[i] = beginTime + i * sweepTime / (POINTS_NUM - 1);
	}




}


/*
* @brief 使用双圆弧拟合算法拟合贝塞尔曲线段
* @param controlPoints	贝塞尔曲线的控制点
* @param beginTime		起点时间(默认为0.0)
* @param endTime		终点时间(默认为1.0)
* @param TOLERENCE		容差，默认为1.0
* @param POINTS_NUM		在拟合的双圆弧上POINTS_NUM个点计算拟合的误差，默认为100
* @param MIN_LENGTH		拟合圆弧的最小长度，默认为1e2，超过这个长度的双圆弧用更多的点计算误差
***/
void BezierToBiArc::fromBezierToBiArc(const std::vector<Vector2d>& controlPoints, const double& beginTime, const double& endTime,const double& TOLERENCE)
{
	assert(endTime > beginTime);

	double splitTime = -1.0;
	const int MAX_ITERATION = 10;
	static const Eigen::Rotation2D<double> antiClockWise(PI / 2);
	static const Eigen::Rotation2D<double> clockWise(-PI / 2);

	Eigen::Vector2d PointA = PointAt(controlPoints, beginTime);
	Eigen::Vector2d PointB = PointAt(controlPoints, endTime);
	Eigen::Vector2d PAPB = (PointB - PointA);
	double dDist = PAPB.norm();

	Eigen::Vector2d leftTangent = computeDerivative(controlPoints, beginTime).normalized();
	Eigen::Vector2d rightTangent = computeDerivative(controlPoints, endTime).normalized();

	//! 两个圆弧的转角总和
	double sweepAngle = std::abs(std::atan2(crossProduct(leftTangent, rightTangent), leftTangent.dot(rightTangent)));
	double Phi = vectorAngle(Eigen::Vector2d(1, 0), leftTangent);
	double Alpha = vectorAngle(leftTangent, PAPB);
	double Beta = vectorAngle(PAPB, rightTangent);



	bool bOppositeSign = std::signbit(Alpha) != std::signbit(Beta);
	//if (bOppositeSign)
	//{
	//	fromBezierToBiArc(controlPoints, beginTime, (beginTime + endTime)/2, TOLERENCE);
	//	fromBezierToBiArc(controlPoints, (beginTime + endTime)/2,endTime , TOLERENCE);
	//	return;
	//}

	double Theta = bOppositeSign ? (3 * Alpha - Beta) / 2 : Alpha;


	//! NOTE: 这里的限制条件是Theta的范围为[-sweepAngle,sweepAngle]
	//! 如果不对自由变量加以限制会怎么样
	//! 
	//Theta = std::clamp(Theta, -sweepAngle, sweepAngle);

	//! 如果自由变量Theta和圆弧的转角sweepAngle相差不大，则调整Theta
	//if (std::abs(Theta) > sweepAngle)
	//{
	//	Theta = Theta > 0 ? Theta - 0.1 : Theta + 0.1;
	//}

	if (std::abs(Alpha + Beta) < 1e-6 || std::abs(Alpha - Beta) < 1e-6)
	{
		std::cout << "ALPHA IS SO CLOSE TO BETA ,ADJUST BETA WITH A DISTURBANCE ." << std::endl;
		Beta += 1e-5;
		//continue;
	}
	for(int i = 0;i<MAX_ITERATION;i++)
	{
		Eigen::Vector2d tangentC(std::cos(Theta + Phi), std::sin(Theta + Phi));
		Eigen::Vector2d normalA = crossProduct(leftTangent, tangentC) > 0 ? antiClockWise * leftTangent : clockWise * leftTangent;
		Eigen::Vector2d normalB = crossProduct(tangentC, rightTangent) > 0 ? antiClockWise * rightTangent : clockWise * rightTangent;
		Eigen::Vector2d normalC = crossProduct(leftTangent, tangentC) > 0 ? antiClockWise * tangentC : clockWise * tangentC;


		double radiusA = dDist * std::sin((Beta - Alpha + Theta) / 2);
		double radiusB = dDist * std::sin((2 * Alpha - Theta) / 2);
		radiusA /= (2 * std::sin((Alpha + Beta) / 2) * std::sin(Theta / 2));
		radiusB /= (2 * std::sin((Alpha + Beta) / 2) * std::sin((Alpha + Beta - Theta) / 2));
		radiusA = std::abs(radiusA);
		radiusB = std::abs(radiusB);


		Eigen::Vector2d PointC = PointA + radiusA * (normalA - normalC);
		Eigen::Vector2d centerA = PointA + radiusA * normalA;
		Eigen::Vector2d centerB = PointB + radiusB * normalB;

		BiArc biArc(PointA, centerA, PointC, centerB, PointB);

		auto errorVec = computeMaxError(biArc, controlPoints, beginTime, endTime, splitTime);
		assert(errorVec.size() == 5);

		//! 如果计算的最大误差小于容差TOLERENCE，则认为拟合成功，返回拟合的双圆弧
		if (*(errorVec.rbegin()) < TOLERENCE)
		{
			m_vecBiArc.emplace_back(biArc);
			return;
		}

		//! 如果计算的最大误差大于容差TOLERENCE，我们希望经过一些优化，
		//! 以便尽量少的使用双圆弧来拟合贝塞尔曲线
		double dErrorA = errorVec[0] + errorVec[1];
		double dErrorB = errorVec[2] + errorVec[3];

		double dSimilarTerms = (4 * std::sin(Alpha + Beta) / 2) / (dDist * std::sin((Alpha - Beta) / 2));
		double deltaThetaA = dErrorA * std::pow(std::sin(Theta / 2), 2) * dSimilarTerms;
		double deltaThetaB = dErrorB * std::pow(std::sin((Alpha + Beta - Theta) / 2), 2) * dSimilarTerms;

		double deltaTheta = std::min(deltaThetaA, -deltaThetaB);
		deltaTheta = std::clamp(deltaTheta, -0.1, 0.1);

		Theta += deltaTheta;
	}

	fromBezierToBiArc(controlPoints,beginTime, splitTime,TOLERENCE);
	fromBezierToBiArc(controlPoints,splitTime, endTime, TOLERENCE);
}

void BezierToBiArc::fromBezierToBiArc_IP(const std::vector<Vector2d>& controlPoints, const double& beginTime, const double& endTime, const double& TOLERENCE)
{
	assert(controlPoints.size() == 4 && endTime > beginTime);

	static const Eigen::Rotation2D<double> antiClockWise(PI / 2);
	static const Eigen::Rotation2D<double> clockWise(-PI / 2);

	//! 计算点A、B的位置、切向、法向
	Eigen::Vector2d PointA = PointAt(controlPoints, beginTime);
	Eigen::Vector2d TangnetA = computeDerivative(controlPoints, beginTime).normalized();
	Eigen::Vector2d NormalA = antiClockWise * TangnetA;

	Eigen::Vector2d PointB = PointAt(controlPoints, endTime);
	Eigen::Vector2d TangentB = computeDerivative(controlPoints, endTime).normalized();
	Eigen::Vector2d NormalB = antiClockWise * TangentB;

	Eigen::Vector2d G = lineIntersection(PointA, TangnetA, PointB, TangentB);
	if (G==Eigen::Vector2d(0,0))
	{
		double midTime = (beginTime + endTime) / 2;
		fromBezierToBiArc_IP(controlPoints, beginTime, midTime, TOLERENCE);
		fromBezierToBiArc_IP(controlPoints, midTime, endTime, TOLERENCE);
		return;
	}
	Eigen::Vector2d PointC = innerPoint(PointA, PointB, G);

	
	Eigen::Vector2d mid_AC = (PointC + PointA) / 2;
	Eigen::Vector2d perPend_AC = antiClockWise * (PointC - PointA).normalized();
	Eigen::Vector2d circleCenterA = lineIntersection(PointA, NormalA, mid_AC, perPend_AC);


	Eigen::Vector2d mid_BC = (PointC + PointB) / 2;
	Eigen::Vector2d perPend_BC = antiClockWise * (PointB - PointC).normalized();
	Eigen::Vector2d circleCenterB = lineIntersection(PointB, NormalB, mid_BC, perPend_BC);


	{
		//! 判断圆心A、圆心B、点C是否共线，如果不共线，计算错误
		Eigen::Vector2d centerAtoC = PointC - circleCenterA;
		Eigen::Vector2d centerBtoC = PointC - circleCenterB;
		if (std::abs(centerAtoC.x() * centerBtoC.y() - centerAtoC.y() * centerBtoC.x()) > 1e-4)
		{
			std::cout << "POINT_C, CENTER_A, CENTER_B ARE NOT COLLINEAR, ERROR."<<std::endl;
		}
	}


	double radiusA = (PointA - circleCenterA).norm();
	double radiusB = (PointB - circleCenterB).norm();

	double radiusAPrime = (PointC - circleCenterA).norm();
	double radiusBPrime = (PointC - circleCenterB).norm();

	{
		//! 判断圆心与圆弧起止点的距离是否相等，如果不相等，计算错误
		if (std::abs(radiusA - radiusAPrime) > 1e-4 || std::abs(radiusB - radiusBPrime) > 1e-4)
		{
			std::cout << "RADIUS_A, RADIUS_B, RADIUS_APRIME, RADIUS_BPRIME ARE NOT EQUAL, ERROR." << std::endl;
		}
	}


	BiArc biArc(PointA, circleCenterA, PointC, circleCenterB, PointB);

	//! 计算双圆弧的最大误差
	//! 计算双圆弧的有向误差时，如果双圆弧的LengthRatio不接近与1或者0
	//! 而且某一个圆弧的累积有向误差比较小，说明某一个圆弧拟合的比较好，或许可以保留
	double maxErrorIndex = -1.0;
	std::vector<double> errorVec = computeMaxError(biArc, controlPoints, beginTime, endTime, maxErrorIndex);
	if (*(errorVec.rbegin()) < TOLERENCE)
	{
		m_vecBiArc.emplace_back(biArc);
	}
	else
	{
		fromBezierToBiArc_IP(controlPoints, beginTime, maxErrorIndex, TOLERENCE);
		fromBezierToBiArc_IP(controlPoints, maxErrorIndex, endTime, TOLERENCE);
	}
	return;
}


/*
* @brief 计算两个向量的夹角,逆时针为正，顺时针为负
* @param v1	向量1
* @param v2	向量2
* @return	夹角的弧度值
***/
inline double BezierToBiArc::vectorAngle(const Vector2d& v1, const Vector2d& v2)
{
	double dotProduct = v1.dot(v2);
	double crossProduct = v1.x() * v2.y() - v1.y() * v2.x();

	return std::atan2(crossProduct, dotProduct);
}

/*
* @brief 计算两个向量的叉乘
* @param v1	向量1
* @param v2	向量2
* @return	叉乘的结果
***/
inline double BezierToBiArc::crossProduct(const Vector2d& v1, const Vector2d& v2)
{
	return v1.x() * v2.y() - v1.y() * v2.x();
}

/*
* @brief 计算双圆弧的最大误差
* @param biarc			拟合的双圆弧
* @param controlPoints	贝塞尔曲线的控制点
* @param beginTime		起点时间
* @param endTime		终点时间
* @param maxErrorIndex	最大误差对应的控制点索引
* @return				误差向量，总共五个元素
*                       errorVec[0] 代表圆弧A的正误差
*                       errorVec[1] 代表圆弧A的负误差
*                       errorVec[2] 代表圆弧B的正误差
*                       errorVec[3] 代表圆弧B的负误差
*                       errorVec[4] 代表双圆弧AB的最大误差
* @note					将双圆弧的弧长纳入误差计算的考量范围内，如果弧长比较长就用100个点，
* 						如果弧长比较短就用适量的点，误差计算的点越多，计算误差越精确。
* @default				默认计算的点数为100个
***/
vector<double> BezierToBiArc::computeMaxError(BiArc& biarc, const vector<Vector2d>& controlPoints, const double& beginTime, const double& endTime, double& maxErrorIndex)
{
	assert(controlPoints.size() == 4);

	vector<double> errorVec(5,0);
	double maxError = 0.0;
	double error = 0.0;
	double totalArcLength = biarc.totalLength;
	const int POINTS_NUM = 100;
	double sweepTime = endTime - beginTime;

	Vector2d Vec = biarc.arc2.centerPoint - biarc.arc1.centerPoint;

	for (size_t i = 0; i < POINTS_NUM; i++)
	{
		double t1 = i * 1.0 / POINTS_NUM;
		Vector2d bezierPoint = PointAt(controlPoints, beginTime + t1 * sweepTime);

		//! 判断贝塞尔曲线上的点在计算误差时属于圆弧A还是圆弧B,依据弧长去判断
		if (t1 <= biarc.lengthRatio)
		{
			error = (bezierPoint - biarc.arc1.centerPoint).norm() - biarc.arc1.radius;
			error > 0 ? errorVec[0] += error : errorVec[1] += error;
		}
		else
		{
			error = (bezierPoint - biarc.arc2.centerPoint).norm() - biarc.arc2.radius;
			error > 0 ? errorVec[2] += error : errorVec[3] += error;
		}
		if (std::abs(error) > maxError)
		{
			maxError = std::abs(error);
			maxErrorIndex = t1;
		}
	}
	errorVec[4] = maxError;

	maxErrorIndex = std::clamp(maxErrorIndex*sweepTime + beginTime, beginTime, endTime);


	//for(size_t i = 0;i<POINTS_NUM;i++)
	//{
	//	double t1 = i * 1.0 / POINTS_NUM;
	//	double t2 = beginTime + t1 * sweepTime;
	//	Vector2d arcPoint = biarc.pointAt(t1);
	//	Vector2d bezierPoint = PointAt(controlPoints, t2);
	//	double error = (arcPoint - bezierPoint).norm();
	//	if (error > maxError)
	//	{
	//		maxError = error;
	//		maxErrorIndex = i;
	//	}
	//}


	return errorVec;
}



/*
* @brief 计算贝塞尔曲线曲率突变的点
* @param controlPoints	贝塞尔曲线的控制点
* @return 	double		曲率突变的点的时间。默认均为2，表示没有曲率变化的点
***/
double BezierToBiArc::getCurvateMutations(const std::vector<Vector2d>& controlPoints)
{
	double t1 = 2.0;
	double t2 = 2.0;

	Eigen::Vector2d A = (-controlPoints[0] + 3 * controlPoints[1] - 3 * controlPoints[2] + controlPoints[3]);
	Eigen::Vector2d B = (controlPoints[0] - 2 * controlPoints[1] + controlPoints[2]);
	Eigen::Vector2d C = (-controlPoints[0] + controlPoints[1]);

	double a1 = A.x();
	double a2 = A.y();
	double b1 = B.x();
	double b2 = B.y();
	double c1 = C.x();
	double c2 = C.y();

	double delta = b1*b1-4*a1*c1;

	//! 如果判别式Δ小于0，则没有实根，曲率没有突变的地方
	//! 如果 a1==0 && b1==0 ，则方程 a1*x^2+b1*x+c1 = 0 变成了c1 = 0，
	//! 当 C1 != 0 时，方程出错，当 C1 == 0 时，整个[0,1]都使得 a1*x^2+b1*x+c1 = 0。
	//! 这个时候一阶导恒为 0，二阶导也恒为 0，从而曲率为 0，也就是说贝塞尔曲线是一条垂直于X轴的直线。
	//! 但是我们已经将直线检测出来了，也就不会有这种情况发生。同样可以类似于方程 a2*x^2+b2*x+c2 = 0 的情况。
	if (delta < 0 || (a1==0 && b1==0)) {
		return 2;
	}

	t1 = a1 == 0 ? -c1 / b1 : (-b1 + sqrt(delta)) / (2 * a1);
	t2 = a1 == 0 ? 2 : (-b1 - sqrt(delta)) / (2 * a1);


	//! 这里还是要判断t1,t2的范围，因为t2 == 2有可能会误传为曲率突变的点
	if (t1 > 0.0 && t1 < 1.0)
	{
		t1 = a2 * t1 * t1 + b2 * t1 + c2 == 0 ? t1 : 2;
	}
	else
	{
		t1 = 2;
	}

	if (t2 > 0.0 && t2 < 1.0)
	{
		t2 = a2 * t2 * t2 + b2 * t2 + c2 == 0 ? t2 : 2;
	}
	else
	{
		t2 = 2;
	}

	//! 判断t1,t2的顺序，输出小的那个
	if (t1 > t2)
	{
		swap(t1, t2);
	}

	return t1;
}

