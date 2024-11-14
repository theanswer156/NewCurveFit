#pragma once
#ifndef BIARC_H
#define BIARC_H
#include <vector>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

static const double PI = 3.14159265358979323846;


/*
* @brief				圆弧类
* @param beginPoint		圆弧起点
* @param centerPoint	圆弧圆心
* @param endPoint		圆弧终点
* @param radius			圆弧半径
* @param startAngle		圆弧起始角度(弧度为单位)
* @param sweepAngle		圆弧掠过角度(弧度为单位),正值表示逆时针画圆弧，负值表示顺时针画圆弧
* @param length			圆弧长度
* @note					统一为逆时针画圆弧，也就是说要求 startAngle < endAngle
*						范围是[0, 2*pi)
*/
class Arc
{
public:
	Arc() {};
	Arc(Eigen::Vector2d beginPoint, Eigen::Vector2d centerPoint, Eigen::Vector2d endPoint);
	~Arc();
	friend std::ostream& operator<<(std::ostream& os, const Arc& arc);
	inline Eigen::Vector2d pointAt(double t);
public:
	double length = 0.0;
	Eigen::Vector2d beginPoint,centerPoint,endPoint;
	double radius=0.0;
	double startAngle = 0.0;
	double endAngle = 0.0;
	double sweepAngle = 0.0;
};

/*
* @brief 双圆弧类, 由两个圆弧构成,两圆弧在公共点处G1连续
* @param arc1 第一个圆弧
* @param arc2 第二个圆弧
* @param lengthRatio 圆弧1与两个圆弧弧长之和的比值，在计算双圆弧上的点时使用
* @param arcBase_1 如果计算的点在圆弧1上，要将传入的t值乘以arcBase_1，得到实际的t值
* @param arcBase_2 如果计算的点在圆弧2上，要将传入的t值减去lengthRatio再乘以arcBase_2，得到实际的t值
* @note 要求 arc1.endPoint == arc2.beginPoint
* @note 默认第一个圆弧在前，第二个圆弧在后。
*
*/
class BiArc
{
public:
	BiArc():totalLength(0.0), lengthRatio(0.0),arcBase_1(0.0), arcBase_2(0.0) {};
	BiArc(const Arc& arc1, const Arc& arc2);
	BiArc(Eigen::Vector2d beginPoint1, Eigen::Vector2d centerPoint1, Eigen::Vector2d endPoint1, Eigen::Vector2d centerPoint2, Eigen::Vector2d endPoint2);
	~BiArc();

	friend std::ostream& operator<<(std::ostream& os, const BiArc& biarc);

	inline Eigen::Vector2d pointAt(const double& t);
public:
	Arc arc1, arc2;
	double totalLength;
	double lengthRatio;
	double arcBase_1, arcBase_2;
};

#endif // !BIARC_H

