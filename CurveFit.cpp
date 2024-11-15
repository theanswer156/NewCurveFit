#include "CurveFit.h"
#include "BezierCurve.h"
#include "BiArc.h"

CurveFit::CurveFit()
{
	getData();
	denseData();
	//this->m_vecSrcData = this->m_vecRawData;
	doComputeTangent_LSE();
	doLineDetect();
	doBiArcFit();
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
			//std::cout << line << std::endl;
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
			if (index != string::npos && num == 8)
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


void CurveFit::doComputeTangent_LSE()
{
	assert(!m_vecSrcData.empty() && m_vecTangentData.empty());
	m_vecTangentData.resize(m_vecSrcData.size(),Eigen::Vector2d());
	double a = 0;
	double b = 0;
	double X = 0;
	double Y = 0;
	double XY = 0;
	double XX = 0;
	int WINDOWSIZE = 9;
	static const int HALFWINDOWSIZE = std::floor(WINDOWSIZE / 2);
	size_t iBegin = 0;
	size_t iEnd = static_cast<size_t> (WINDOWSIZE - 1);
	size_t iMid = static_cast<size_t> (std::floor(iBegin + iEnd) / 2);

	for (size_t index = 0; index < iMid; index++)
	{
		m_vecTangentData[index] = (m_vecSrcData[index + 1] - m_vecSrcData[index]).normalized();
	}

	for (size_t index = 0; index < WINDOWSIZE; index++)
	{
		X += m_vecSrcData[index].x();
		Y += m_vecSrcData[index].y();
		XY += m_vecSrcData[index].x() * m_vecSrcData[index].y();
		XX += m_vecSrcData[index].x() * m_vecSrcData[index].x();
	}

	while (iEnd < m_vecSrcData.size())
	{
		bool bAbNormal = false;
		if (iBegin)
		{
			X -= (m_vecSrcData[iBegin - 1].x());
			X += (m_vecSrcData[iEnd].x());
			Y -= (m_vecSrcData[iBegin - 1].y());
			Y += (m_vecSrcData[iEnd].y());
			XX -= (m_vecSrcData[iBegin - 1].x() * m_vecSrcData[iBegin - 1].x());
			XX += (m_vecSrcData[iEnd].x() * m_vecSrcData[iEnd].x());
			XY -= (m_vecSrcData[iBegin - 1].x() * m_vecSrcData[iBegin - 1].y());
			XY += (m_vecSrcData[iEnd].x() * m_vecSrcData[iEnd].y());
		}


		if (std::abs(WINDOWSIZE * XX - X * X) < 1e-4 || std::abs(WINDOWSIZE * XY - X * Y) < 1e-4)
		{
			//std::cout << "IMID : " << iMid << "\t" << "IBEGIN : " << iBegin << "\t" << "IEND : " << iEnd << std::endl;
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
			Eigen::Vector2d lineVector = m_vecSrcData[iEnd] - m_vecSrcData[iMid];
			if (bAbNormal)
			{
				bool bXDirection = std::abs(lineVector.x()) > std::abs(lineVector.y());
				double dValue = 0.0;
				if (bXDirection)
				{
					dValue = lineVector.x() > 0 ? 1 : -1;
					m_vecTangentData[iMid] = Eigen::Vector2d{ dValue,0 };
				}
				else
				{
					dValue = lineVector.y() > 0 ? 1 : -1;
					m_vecTangentData[iMid] = Eigen::Vector2d{ 0, dValue };
				}
			}
			else
			{
				//!	默认Y的分量为正
				double modulus = std::sqrt(1 + a * a);
				bool bSameDirect = lineVector.x() * a + lineVector.y() >= 0;
				//m_vecTangentData[iMid] = bSameDirect ? Point{ 1 / modulus,a / modulus } : Point{ -1 / modulus,-a / modulus };
				m_vecTangentData[iMid] = Eigen::Vector2d{ 1 / modulus,a / modulus };
				double angle = std::atan2(m_vecTangentData[iMid].x() * lineVector.y() - m_vecTangentData[iMid].y() * lineVector.x(), m_vecTangentData[iMid].dot(lineVector));
				m_vecTangentData[iMid] = std::abs(angle) <= PI / 2 ? m_vecTangentData[iMid] : -m_vecTangentData[iMid];

			}
		}
		++iBegin;
		++iMid;
		++iEnd;
	}
}


/*
* @brief 判断曲线的起始点和终止点的切向是否需要改变
* @param begin 起始点索引
* @param end 终止点索引
* @note 起始点就看后面的点的切向与他的切向变化关系，如果变化大，使用后向差分公式计算切向
*       终止点就看前面的点的切向与他的切向变化关系，如果变化大，使用前向差分公式计算切向
***/
void CurveFit::isTangentChange(const int& begin, const int& end)
{
	int WINDOWSIZE = 9;
	static const int HALFWINDOWSIZE = std::floor(WINDOWSIZE / 2);
	const double ANGEL_THRESHOLD = PI/8;
	//! 判断起始点是否需要使用后向差分公式计算切向
	Eigen::Vector2d tangentBegin = this->m_vecTangentData[begin];
	Eigen::Vector2d tangent_1 = this->m_vecTangentData[begin + HALFWINDOWSIZE];
	double sweepAngle_Begin = std::atan2(tangent_1.x() * tangentBegin.y() - tangent_1.y() * tangentBegin.x(), tangent_1.dot(tangentBegin));

	if (begin - HALFWINDOWSIZE > 0)
	{
		Eigen::Vector2d tangent_2 = this->m_vecTangentData[begin - HALFWINDOWSIZE];
		double sweepAngle = std::atan2(tangent_2.x() * tangentBegin.y() - tangent_2.y() * tangentBegin.x(), tangent_2.dot(tangentBegin));
		sweepAngle_Begin+=sweepAngle;
	}
	if (sweepAngle_Begin > ANGEL_THRESHOLD)
	{
		this->m_vecTangentData[begin] = (this->m_vecSrcData[begin + HALFWINDOWSIZE] - this->m_vecSrcData[begin]).normalized();
	}



	Eigen::Vector2d tangentEnd = this->m_vecTangentData[end];
	Eigen::Vector2d tangent_3 = this->m_vecTangentData[end - HALFWINDOWSIZE];
	double sweepAngle_End = std::atan2(tangent_3.x() * tangentEnd.y() - tangent_3.y() * tangentEnd.x(), tangent_3.dot(tangentEnd));

	if (end + HALFWINDOWSIZE < m_vecSrcData.size())
	{
		Eigen::Vector2d tangent_4 = this->m_vecTangentData[end + HALFWINDOWSIZE];
		double sweepAngle = std::atan2(tangent_4.x() * tangentEnd.y() - tangent_4.y() * tangentEnd.x(), tangent_4.dot(tangentEnd));
		sweepAngle_End += sweepAngle;
	}

	if (sweepAngle_End > ANGEL_THRESHOLD)
	{
		this->m_vecTangentData[end] = (this->m_vecSrcData[end] - this->m_vecSrcData[end - HALFWINDOWSIZE]).normalized();
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
	const int LENGTH_THRESHOLD = 6;
	size_t srcDataSize = m_vecSrcData.size();
	m_bvecStraightFlags.resize(m_vecSrcData.size(),false);
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
		//! 判断第一个曲线段是直线还是曲线，如果vecCurveCutIndex[2]!= 0，说明第一个曲线段是曲线
		//! 因此我们要先把第一个曲线段的直线部分放入m_vecLineSegments中，
		//! 然后就能按照下面的方式依次处理直线和曲线交替出现的情况
		if (vecCurveCutIndex[Index + 1])
		{
			++Index;
			size_t iCurveEnd = vecCurveCutIndex[Index];
			this->m_vecCurveSegments.emplace_back(std::make_pair(0, iCurveEnd));
		}
		while (Index < vecCurveCutIndex.size() - 2)
		{
			size_t iStraightBegin = vecCurveCutIndex[Index];
			++Index;
			while (Index < vecCurveCutIndex.size() - 1 && vecCurveCutIndex[Index] == 0)
			{
				++Index;
			}

			size_t iStraightEnd = vecCurveCutIndex[Index];
			//! 先前假设是先直线然后紧接着是曲线，但是可能是先曲线后直线，
			//! 所以在前面最好要判断首先出现的是曲线还是直线，然后再进行处理
			if (Index < vecCurveCutIndex.size() )
			{
				this->m_vecLineSegments.emplace_back(std::make_pair(iStraightBegin, iStraightEnd));
			}

			++Index;

			if (Index < vecCurveCutIndex.size() )
			{
				size_t iCurveBegin = iStraightEnd;
				size_t iCurveEnd = vecCurveCutIndex[Index];
				this->m_vecCurveSegments.emplace_back(std::make_pair(iCurveBegin, iCurveEnd));
			}
		}
	}


}

void CurveFit::doBiArcFit()
{
	for (const auto& pair : m_vecCurveSegments)
	{
		int iBeginIndex = pair.first;
		int iEndIndex = pair.second;

		isTangentChange(iBeginIndex, iEndIndex);

		std::vector<Eigen::Vector2d> subVec(m_vecSrcData.begin() + iBeginIndex, m_vecSrcData.begin() + iEndIndex + 1);
		std::vector<Eigen::Vector2d> subTangnent(m_vecTangentData.begin() + iBeginIndex, m_vecTangentData.begin() + iEndIndex + 1);

		//! 考虑到曲线起始点和终止点切向对贝塞尔曲线拟合的影响很大，
		//! 我们希望计算出起始点和终止点的一个邻域内的切线的变化计算出来，
		//! 只看变化的方向，不看变化的大小，两两相比较，
		//! 如果变化的很大说明采用LSE的方法计算出来的切向有问题，
		//! 从而我们选取前向或者后向差分公式来计算切向
		m_vecCurves.emplace_back(subVec);
		m_vecTangent.emplace_back(subTangnent);
	}
}
