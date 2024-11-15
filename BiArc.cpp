#include "BiArc.h"
Arc::Arc(Eigen::Vector2d beginPoint, Eigen::Vector2d centerPoint, Eigen::Vector2d endPoint) :
	beginPoint(beginPoint), centerPoint(centerPoint), endPoint(endPoint)
{
	//! 断言beginPoint到centerPoint和endPoint到centerPoint的距离几乎没有什么误差
	double distance_1 = (beginPoint - centerPoint).norm();
	double distance_2 = (endPoint - centerPoint).norm();
	
	double distance_error = std::abs(distance_1 - distance_2);
	std::cout << " ARC distance_error: " << distance_error << std::endl;

	//assert(std::abs((beginPoint - centerPoint).norm() - (endPoint - centerPoint).norm()) < 1e-4);
	radius = (centerPoint - beginPoint).norm();
	Eigen::Vector2d centerToBegin = beginPoint - centerPoint;
	Eigen::Vector2d centerToEnd = endPoint - centerPoint;
	startAngle = atan2(centerToBegin.y(), centerToBegin.x());
	endAngle = atan2(centerToEnd.y(), centerToEnd.x());
	if(startAngle < 0) startAngle += 2 * PI;
	if (endAngle < 0) endAngle += 2 * PI;
	if (startAngle > endAngle)
	{
		std::swap(startAngle, endAngle);
	}
	sweepAngle = endAngle - startAngle;
	length = std::abs(sweepAngle * radius);
	assert(sweepAngle >= 0);
}
Arc::~Arc()
{
}

inline Eigen::Vector2d Arc::pointAt(double t)
{
	double angle = startAngle + t * sweepAngle;
 
	return centerPoint + radius * Eigen::Vector2d(cos(angle), sin(angle));
}





BiArc::BiArc(const Arc& arc1,const Arc& arc2):arc1(arc1), arc2(arc2)
{
	totalLength = arc1.length + arc2.length;
	lengthRatio = arc1.length / totalLength;
	arcBase_1 = 1 / lengthRatio;
	arcBase_2 = 1 / (1 - lengthRatio);
}

/*
* @brief: 五点构造一个双圆弧类
* @param: beginPoint1 第一条圆弧的起点
* @param: centerPoint1 第一条圆弧的圆心
* @param: endPoint1 第一条圆弧的终点
* @param: centerPoint2 第二条圆弧的圆心
* @param: endPoint2 第二条圆弧的终点
***/
BiArc::BiArc(Eigen::Vector2d beginPoint1, Eigen::Vector2d centerPoint1, Eigen::Vector2d endPoint1, Eigen::Vector2d centerPoint2, Eigen::Vector2d endPoint2)
{
	arc1 = Arc(beginPoint1, centerPoint1, endPoint1);
	arc2 = Arc(endPoint1, centerPoint2, endPoint2);
	totalLength = arc1.length + arc2.length;
	lengthRatio = arc1.length / totalLength;
	arcBase_1 = 1 / lengthRatio;
	arcBase_2 = 1 / (1 - lengthRatio);
}

BiArc::~BiArc()
{
}
/*
* @brief:	计算双圆弧上某一点的坐标
* @param: t 双圆弧上的比例参数，范围为[0,1]
* @return:	返回双圆弧上t比例处的点坐标
* @note:	如果t>lengthRatio,返回的是第一条圆弧上的点，否则返回的是第二条圆弧上的点
*/
inline Eigen::Vector2d BiArc::pointAt(const double& t)
{
	return t <= lengthRatio? arc1.pointAt(arcBase_1 * t) : arc2.pointAt((t - lengthRatio) * arcBase_2);
}


/*
*@brief		重载圆弧类的输出流运算符
* @param	os 输出流对象
* @param	arc 圆弧对象
* @return	返回输出流对象
***/
std::ostream& operator<<(std::ostream& os, const Arc& arc)
{
	os << arc.centerPoint.x() << "\t" << arc.centerPoint.y() << "\t" << arc.radius << "\t" << arc.startAngle << "\t" << arc.endAngle << "\t" << arc.sweepAngle << std::endl;
	return os;
}



/*
*@brief		重载双圆弧类的输出流运算符
* @param	os 输出流对象
* @param	arc 圆弧对象
* @return	返回输出流对象
***/
std::ostream& operator<<(std::ostream& os, const BiArc& biarc)
{
	os << biarc.arc1;
	os << biarc.arc2;
	return os;
}
