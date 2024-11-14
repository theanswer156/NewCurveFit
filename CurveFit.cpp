#include "CurveFit.h"
#include "BezierCurve.h"
#include "BiArc.h"

CurveFit::CurveFit()
{
}

CurveFit::~CurveFit()
{
}
/**
 * @brief 从plt文件中读取原始数据,存放到m_vecRawData中
 */
void CurveFit::getData()
{
	size_t size = 100;
	std::vector<std::vector<std::string>> vecString(size, std::vector<std::string>());
	int index = 0;
	std::string pltFilePath = "C:/Users/Zhushengb/source/repos/CurveFit/alphabet.plt";
	std::ifstream pltFile(pltFilePath);

	bool bBeginRead = false;
	bool bEndRead = false;

	if (!pltFile.is_open()) {
		std::cerr << "Error opening file: " << pltFilePath << std::endl;
		return;
	}

	std::string line;
	while (std::getline(pltFile, line)) {
		// 打印文件的每一行
		std::cout << line << "\t" << std::endl;
		if (line.substr(0, 2) == "PU") {
			bBeginRead = true;
		}
		if (bBeginRead && !bEndRead)
		{
			vecString[index].emplace_back(line);
		}
		if ((line.substr(0, 2) == "SP") && bBeginRead)
		{
			std::cout << line << std::endl;
			//std::cout << std::endl;
			//std::cout << std::endl;
			//std::cout << std::endl;
			bBeginRead = false;
			bEndRead = false;
			index++;
		}
		if (index == size - 1)
		{
			break;
		}
		//std::cout << line << std::endl;
	}
	pltFile.close();

	std::cout << "begin output line massage. " << std::endl;
	int num = 0;
	for (const auto& strs : vecString)
	{
		for (const auto& str : strs)
		{
			size_t index = str.find(" ");
			if (index != string::npos && num == 10)
			{
				Eigen::Vector2d vec;
				std::string x = str.substr(2, index - 2);
				std::string y = str.substr(index + 1, str.size() - index - 1);
				vec<<stod(x), stod(y);
				m_vecRawData.emplace_back(vec);
				//std::cout << "X : " << x << "\t" << "Y : " << y << std::endl;
			}
			//std::cout << str << std::endl;
		}
		num++;
		//std::cout << std::endl;
		//std::cout << std::endl;
		//std::cout << std::endl;

	}
}
/*
*	@brief 将原始数据稠密化，存放到m_vecSrcData中
*   @param interval 采样间隔
* 	@note 这里的稠密化是指将原始数据按照一定间隔进行采样，
*		 并将采样得到的点存放到m_vecSrcData中
***/
void CurveFit::denseData()
{
	assert(m_vecRawData.size()&& !m_vecSrcData.size());
	size_t n = this->m_vecRawData.size();
	int iBegin = 0;
	int iEnd = 1;
	const double MIN_DISTANCE = 2.0;
	while (iEnd < n)
	{
		double dist = (this->m_vecRawData[iEnd] - this->m_vecRawData[iBegin]).norm();
		if (dist != 0.0)
		{
			int INSERTNUMBER = dist < MIN_DISTANCE ? 0 : static_cast<int>(std::ceil(dist / MIN_DISTANCE));;
			this->m_vecSrcData.emplace_back(this->m_vecRawData[iBegin]);
			Eigen::Vector2d direc = (this->m_vecRawData[iEnd] - this->m_vecRawData[iBegin]);
			for (int i = 1; i < INSERTNUMBER; i++)
			{
				Eigen::Vector2d insertPoint = this->m_vecRawData[iBegin] + (i * 1.0 / INSERTNUMBER) * direc;
				this->m_vecSrcData.emplace_back(insertPoint);
			}
		}
		++iEnd;
		++iBegin;
	}
}
/**
 * @brief 计算向量叉乘
 * @param vec1 向量1
 * @param vec2 向量2
 * @return 向量叉乘结果
 */
double CurveFit::crossProduct(const Eigen::Vector2d& vec1, const Eigen::Vector2d& vec2)
{
	return vec1.x() * vec2.y() - vec1.y() * vec2.x();
}
radian CurveFit::VectorAngle(const Eigen::Vector2d& vec1, const Eigen::Vector2d& vec2)
{
	radian ans = 0.0;
	double cosAlpha = (vec1.x() * vec2.x() + vec1.y() * vec2.y()) / (vec1.norm() * vec2.norm());
	double sinAlpha = crossProduct(vec1, vec2);
	if (std::abs(sinAlpha) > 1e-8)
	{
		return !std::signbit(cosAlpha) ? std::acos(cosAlpha) : -std::acos(cosAlpha);
	}
	else
	{
		return !std::signbit(cosAlpha) ? 0.0 : PI;
	}
	return ans;
}
/*
*	@brief 对稠密后的原始数据进行直线检测通过计算相邻三个点之间的角度，
*          判断是否为直线，进行直线检测，并将直线、曲线部分划分
*   @param ANGEL_THRESHOLD 直线检测的角度阈值
*   @param LENGTH_THRESHOLD 直线检测的长度阈值
* 
*/
void CurveFit::doLineDetect()
{
	assert(m_vecSrcData.size());
	const radian ANGEL_THRESHOLD = 1e-8;
	const int LENGTH_THRESHOLD = 5;
	size_t srcDataSize = m_vecSrcData.size();
	m_bvecStraightFlags.resize(m_vecSrcData.size());
	//! 遍历m_vecsrcData中的每一个点，判断是否为直线
	for (size_t i = 1; i <= m_vecSrcData.size() - 2; i++)
	{
		Eigen::Vector2d vecBefore = (m_vecSrcData[i - 1] - m_vecSrcData[i]).normalized();
		Eigen::Vector2d vecAfter = (m_vecSrcData[i + 1] - m_vecSrcData[i]).normalized();
		double crossProductResult = crossProduct(vecBefore, vecAfter);
		if (std::abs(crossProductResult) < ANGEL_THRESHOLD)
		{
			//! 如果叉乘的绝对值接近零，说明这相邻的三个点近似在一条直线上
			m_bvecStraightFlags[i] = true;
		}
		else
		{
			m_bvecStraightFlags[i] = false;
		}
	}

	// 在对数据进行初步的直线检测之后，我们还要对数据进行处理，
	// 把连续直线之间出现的小段曲线设置为直线
	// 把连续曲线之间的小段直线设置为曲线 
	// 曲线和直线的最小长度都设置一个阈值




	//! 将直线段的起点和终点标记为直线段
	//! 遍历m_bvecStraightFlags中的每一个点，如果该点为true，则该点为直线段上的点
	//! 先放两个0，最后放上两个m_vecSrcData.size() - 1，这样比较好处理
	vector<size_t> vecCurveCutIndex(2,0);
	size_t iTrueStart = -1;
	int iTrueLength = 0;

	//! 将连续出现的true的情形(直线部分)划分为直线段
	for (size_t i = 0; i < m_bvecStraightFlags.size(); i++)
	{
		if (m_bvecStraightFlags[i])
		{
			if (iTrueStart == -1)
			{
				iTrueStart = i;
			}
			++iTrueLength;
		}
		else
		{
			if (iTrueLength >= LENGTH_THRESHOLD)
			{
				vecCurveCutIndex.emplace_back(iTrueStart);
				vecCurveCutIndex.emplace_back(i);

			}
			iTrueStart = -1;
			iTrueLength = 0;
		}
	}
	//! 处理最后连续出现true的情况
	if (iTrueLength >= LENGTH_THRESHOLD)
	{
		vecCurveCutIndex.emplace_back(iTrueStart);
		vecCurveCutIndex.emplace_back(srcDataSize - 1);
	}

	vecCurveCutIndex.emplace_back(srcDataSize - 1);
	vecCurveCutIndex.emplace_back(srcDataSize - 1);

	std::cout << string(50, '-') << std::endl;
	std::cout << "vecCurveCutIndex size : " << vecCurveCutIndex.size() << std::endl;
	std::cout << string(50, '-') << std::endl;


	//! 在检测出来的直线之间就是我们要找的曲线
	//! 在vector<size_t> vecCurveCutIndex中存放了曲线的分割点的索引
	//! 其中为了方便，我们在vecCurveCutIndex的首尾各放上两个点，这样比较好处理
	int MIN_CURVE_LENGTH = 5;
	//! vector是动态的，所以vecCurveCutIndex.size() - 2也是会变的，因此不能直接用
	for (size_t i = 1; i < vecCurveCutIndex.size() - 2; i += 2)
	{
		if (vecCurveCutIndex[i + 1] - vecCurveCutIndex[i] < MIN_CURVE_LENGTH)
		{
			vecCurveCutIndex[i] = 0;
			vecCurveCutIndex[i + 1] = 0;
		}
	}


	{
		size_t Index = 1;
		while (Index < vecCurveCutIndex.size() - 2)
		{
			size_t iStraightBegin = vecCurveCutIndex[Index];
			++Index;
			while (Index < vecCurveCutIndex.size() - 2 && vecCurveCutIndex[Index] == 0)
			{
				++Index;
			}
			size_t iStraightEnd = vecCurveCutIndex[Index];
			this->m_vecLineSegments.emplace_back(std::make_pair(iStraightBegin, iStraightEnd));

			++Index;
			size_t iCurveBegin = iStraightEnd;
			size_t iCurveEnd = vecCurveCutIndex[Index];
			this->m_vecCurveSegments.emplace_back(std::make_pair(iCurveBegin, iCurveEnd));
		}
	}


}

void CurveFit::doBiArcFit()
{
	assert(m_vecLineSegments.size());
	for (const auto& pair : m_vecCurveSegments)
	{
		int iBeginIndex = pair.first;
		int iEndIndex = pair.second;
		std::vector<Eigen::Vector2d> subVec(m_vecSrcData.begin() + iBeginIndex, m_vecSrcData.begin() + iEndIndex + 1);
		m_vecCurves.emplace_back(subVec);
	}
}
