#pragma once

#ifndef CURVEFIT_H
#define CURVEFIT_H
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cassert>
#include <ctime>
#include <chrono>

//const double PI = 3.14159265358979323846;
using radian = double;
using namespace std;
using namespace Eigen;
/*
* @brief 存储曲线段和直线段的起始点和终止点
* @param  begin: 起始点索引(int)
* @param  end: 终止点索引(int)
* @param  isStraightLine: 是否为直线段(bool)
*/
struct Node
{
    int begin;
    int end;
    bool isStraightLine;
};
struct LineSegment
{
    LineSegment() {}
    LineSegment(const Eigen::Vector2d& beginPoint, const Eigen::Vector2d& endPoint) : beginPoint(beginPoint), endPoint(endPoint) {}
    Eigen::Vector2d beginPoint;
    Eigen::Vector2d endPoint;
};

class CurveFit
{
public:
    CurveFit();
    ~CurveFit();
    void getData();
    void readTxtData();
    void denseData();

    void doComputeTangent_LSE();
    void isTangentChange(const int& begin, const int& end);
    double crossProduct(const Eigen::Vector2d& vec1, const Eigen::Vector2d& vec2);
    radian VectorAngle(const Eigen::Vector2d& vec1, const Eigen::Vector2d& vec2);
    void doLineDetect();
    void doBiArcFit();

    const vector<vector<Eigen::Vector2d>>& getCurves() const { return m_vecCurves; }
    const vector<vector<Eigen::Vector2d>>& getTangent() const { return m_vecTangent; }

private:

    std::vector<Eigen::Vector2d> m_vecSrcData;
    std::vector<Eigen::Vector2d> m_vecTangentData;
    std::vector<Eigen::Vector2d> m_vecRawData;
    std::vector<bool> m_bvecStraightFlags;
    std::vector<pair<int, int>> m_vecLineSegments;
    std::vector<pair<int, int>> m_vecCurveSegments;

    vector<vector<Eigen::Vector2d>> m_vecCurves;
    vector<vector<Eigen::Vector2d>> m_vecTangent;
};
#endif // !CURVEFIT_H
