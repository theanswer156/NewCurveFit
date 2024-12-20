#include "CubicBezier.h"

CubicBezier::CubicBezier(const std::vector<Vector2d>& controlPoints)
{
	BezierCurve::setControlPoints(controlPoints);
}

CubicBezier::CubicBezier(const std::vector<Vector2d>& contralPoints, const std::pair<double, double>& range)
{
	BezierCurve::setControlPoints(contralPoints);
	this->m_vecRange = range;
}

/*
*@brief: 使用ranges库构造CubicBezier对象
*@param: const std::ranges::subrange<vector<Vector2d>::iterator> controlPoints: 控制点的范围
*@param: const std::pair<double, double>& range: 范围
***/
CubicBezier::CubicBezier(const std::ranges::subrange<vector<Vector2d>::const_iterator> controlPoints, const std::pair<double, double>& range)
{
	assert(std::distance(controlPoints.begin(), controlPoints.end()) == 4);

	std::vector<Vector2d> controlPointsVec(controlPoints.begin(), controlPoints.end());
	BezierCurve::setControlPoints(controlPointsVec);
	//for (auto const& contrlPoint : controlPoints)
	//{
	//	
	//}
	this->m_vecRange = range;
}


/*
*@brief: 从DXF文件中将三次贝塞尔曲线分离出来
*@param: const std::string& filename: DXF文件名
*@return: vector<CubicBezier>: 分离出的三次贝塞尔曲线
****/
const vector<CubicBezier> CubicBezier::splitDXF(const std::string& file_path)
{
	const string postfix = ".dxf";
	static const string STARTLOG = "AcDbSpline";
	vector<CubicBezier> ans;

	std::ifstream infile(file_path);
	if (file_path.substr(file_path.size() - postfix.size()) != postfix)
	{
		std::cerr << "Error: file name should end with \" .dxf \". " << std::endl;
		return vector<CubicBezier>();
	}
	if (!infile.is_open())
	{
		std::cerr << "Error: cannot open file " << file_path << std::endl;
	}
	std::string line;
	while (std::getline(infile, line))
	{
		auto pos = line.find(STARTLOG)!= std::string::npos;
		if (pos)
		{
			std::cout << "FOUND STRATLOG, START TO TRANSFORM SPLINE TO BEZIER" << std::endl;
			vector<CubicBezier> temp = transSplineToBezier(infile);
			ans.insert(ans.end(), temp.begin(), temp.end());
		}
	}


	return ans;
}

/*
*@brief: 将三次样条曲线转换为贝塞尔曲线
*@param: const std::ifstream& infile: 三次样条曲线数据文件
*@return: vector<CubicBezier>&: 转换后的三阶贝塞尔曲线
***/
vector<CubicBezier> CubicBezier::transSplineToBezier(std::ifstream& infile)
{
	vector<CubicBezier> ans{};
	static const string ENDLOG = "SPLINE";
	int DEGREE = 0;
	int KNOTSNUM = 0;
	int CONTROLNUM = 0;
	int FITTINGPOINTSNUM = 0;
	//! 读取到AcDbSpline字符串以后就开始读取样条曲线的各类数据了
	static auto ignore_N_line = [&infile](const int n)->void
		{
			for (int i = 0; i < n && infile.good(); ++i)
			{
				infile.ignore(50, '\n');
			}
	};
	//! 读取组码71、72、73所表示的样条曲线的次数、节点向量的长度、控制点的数量
	static auto get_param = [&]() -> void
		{
			infile.ignore(50, '\n');
			string line;
			std::getline(infile, line);
			DEGREE = std::stoi(line);// 组码71
			infile.ignore(50, '\n');
			std::getline(infile, line);
			KNOTSNUM = std::stoi(line);// 组码72
			infile.ignore(50, '\n');
			std::getline(infile, line);
			CONTROLNUM = std::stoi(line);// 组码73
			infile.ignore(50, '\n');
			std::getline(infile, line);
			FITTINGPOINTSNUM = std::stoi(line);// 组码74
	};
	
	//! 跳过前面8行有关组码210、220、230 、70的无用数据
	ignore_N_line(8); 
	//! 读取组码71、72、73所表示的样条曲线的次数、节点向量的长度、控制点的数量
	//! 不过这里是默认DXF遵循的规范是AC1032，如果是AC1015，组码72为0
	get_param();
	assert(DEGREE == 3 && "Error: the degree of spline curve should be 3.");

	//! 如果DXF文件遵循的文件规范是AC1015，则KNOTSNUM可能为0，
	//! 此时需要根据DEGREE和控制点的数量计算KNOTSNUM
	if (!KNOTSNUM)
	{
		KNOTSNUM = DEGREE + KNOTSNUM + 1;
	}
	vector<Vector2d> controlPoints(CONTROLNUM,Eigen::Vector2d::Zero());
	//! 注意这里记录的是消重之后的节点向量   所以容量是CONTROLNUM/3
	//! 使用标准库中的allocator有错误  显示没有construct成员
	vector<double> knots(KNOTSNUM, 0.0);

	//! 跳过后面3行有关组码74、42、43的无用数据
	if (FITTINGPOINTSNUM)
	{
		ignore_N_line(6);
	}
	else
	{
		ignore_N_line(4);
	}

	//! 依据组码72、73的值读取节点向量、控制点坐标
	//! 
	//! 
	//! 如果这里的组码74不为0，说明含有拟合点，后面的组码43就会存在然后扰乱整个方法
	//! 所以这里需要判断一下组码74是否为0，如果不为0，则跳过后面的组码43
	//! 
	{
		string line;
		//! 读取节点向量
		for (int i = 0; i < knots.size(); ++i)
		{
			ignore_N_line(1);
			std::getline(infile, line);
			knots[i] = std::stod(line);
		}

		//! 开始读取控制点坐标
		for (int i = 0; i < CONTROLNUM; ++i)
		{
			ignore_N_line(1);
			std::getline(infile, line);
			double x = std::stod(line);
			ignore_N_line(1);
			std::getline(infile, line);
			double y = std::stod(line);
			controlPoints[i] = Vector2d(x, y);
			ignore_N_line(2);
		}
	}


	{
		std::ofstream outfile("C:\\Users\\Zhushengb\\Desktop\\vector_graphic\\knots_008.txt", std::ios::app | std::ios::out);
		if (outfile.is_open())
		{
			for (size_t i = 1; i < knots.size(); i++)
			{
				outfile << "TIME DIFFERENCRE: " << i << "\t" << knots[i] << "\t" << knots[i - 1] << "\t is \t" << knots[i] - knots[i - 1] << std::endl;
			}
			outfile << std::endl;
			outfile << std::endl;
			outfile << std::endl;
			outfile << std::endl;
			outfile.close();
		}
		else
		{
			std::cerr << "Error: cannot open file" << std::endl;
		}


	}


	//! 将读取到的节点向量、控制点坐标转换为三次贝塞尔曲线
	//! Boehm算法一次只能插入一个节点，但是可以算出多个控制点，可是算出来的控制点要插入到控制点向量中
	//! 这就涉及到对vector的频繁插入，如果要插入的节点非常多，效率会很低
	//! 所以这里可以选择奥斯陆算法，一次插入N个节点，但是每次都算一个控制点，这样我们可以将原有的控制点vector和计算出来的vector
	//! 还是写一个函数判断样条函数是分段贝塞尔形式还是端点插值的非均匀三次B样条曲线形式做不同的处理
	//! 
	//! TODO: 后续将两种处理能不能做成一样的处理方式
	bool isUniform = CubicBezier::isUniform(knots);
	if (isUniform)
	{
		uniformSplineToBezier(knots, controlPoints, ans);
	}
	else
	{
		NonuniformSplineToBezier(knots, controlPoints, ans);
	}

	{
		std::ofstream outfile("C:\\Users\\Zhushengb\\source\\repos\\newCurveFit\\BezierCurve.txt", std::ios::app | std::ios::out);
		if (outfile.is_open())
		{
			for (auto& bezier : ans)
			{
				outfile << bezier;
			}
			outfile.close();
		}
		else
		{
			std::cerr << "Error: cannot open The BezierCurve.txt file"<< std::endl;
		}
	}

	return ans;
}

/*
*@brief: 将节点向量转换为map,帮助我们判断这样一条曲线是分段贝塞尔形式还是端点插值的非均匀三次B样条曲线形式
* @param: const std::vector<double>& knotVec: 节点向量
* @param: pair<double,int>: <节点值,节点重复度>
* @return: unordered_map<double, int>: 节点向量的map
* @notice: 由于首节点和尾节点的重复度都是4，所以我们不需要记录
***/
bool CubicBezier::isUniform(const std::vector<double>& knotVec)
{
	bool ans = true;
	unordered_map<double, int> knot_map;
	for (int i = 4; i < knotVec.size() - 4; i++)
	{
		knot_map[knotVec[i]]++;
	}
	for (auto& knot_pair : knot_map)
	{
		if (knot_pair.second == 1)
		{
			ans = false;
			break;
		}
	}
	return ans;
}

/*
*@brief: 将分段三次贝塞尔形式的三次样条曲线转换为贝塞尔曲线
* @param: const std::vector<double>& knotVec: 节点向量
* @param: const std::vector<Vector2d>& controlPoints: 控制点向量
* @param: std::vector<CubicBezier>& bezierVec: 转换后的贝塞尔曲线向量
***/
void CubicBezier::uniformSplineToBezier(const std::vector<double>& knotVec, const std::vector<Vector2d>& controlPoints, std::vector<CubicBezier>& bezierVec)
{
	for (int i = 3; i < knotVec.size() - 4; i += 3)
	{
		pair<double, double> range(knotVec[i], knotVec[i +1]);
		std::ranges::subrange subrange(controlPoints.begin() + i - 3, controlPoints.begin() + i + 1);
		CubicBezier bezier(subrange, range);
		bezierVec.push_back(bezier);
	}
}


/*
*@brief: 将非均匀三次贝塞尔形式的三次样条曲线转换为贝塞尔曲线，
*		 这个函数只适用于中间节点全是一重的情况，如果要修改成能够同时处理三重节点和一重节点混合的情况，
*		 必须要确定有效的节点区间，即要求[t_j,t_{j+1})，t_{j}<t_{j+1},不能是等于的情况。
*		 这个时候最好使用while遍历，在找到相应的有效节点区间时要判断t_{j},t_{j+1}是否是一重的，
*		 如果是不全是一重的，通过下面的计算得到控制点，如果都是三重节点，直接从B样条的控制点中分离出来即可。
* @param: const std::vector<double>& knotVec: 节点向量
* @param: const std::vector<Vector2d>& controlPoints: 控制点向量
* @param: std::vector<CubicBezier>& bezierVec: 转换后的贝塞尔曲线向量
****/
void CubicBezier::NonuniformSplineToBezier(const std::vector<double>& knotVec, const std::vector<Vector2d>& controlPoints, std::vector<CubicBezier>& bezierVec)
{
	for (size_t j = 3; j < knotVec.size() - 4; ++j)
	{
		pair<double, double> range(knotVec[j], knotVec[j + 1]);
		////! P1计算没有问题
		//assert(knotVec[j + 2] > knotVec[j - 1] && "Error: knot vector is not in ascending order.");
		//Vector2d p1 = ((knotVec[j] - knotVec[j - 1]) * controlPoints[j - 1] + (knotVec[j + 2] - knotVec[j]) * controlPoints[j - 2]) / (knotVec[j + 2] - knotVec[j - 1]);

		//assert(knotVec[j + 3] > knotVec[j] && "Error: knot vector is not in ascending order.");
		//Vector2d p2 = ((knotVec[j + 1] - knotVec[j]) * controlPoints[j] + (knotVec[j + 3] - knotVec[j + 1]) * controlPoints[j - 1]) / (knotVec[j + 3] - knotVec[j]);

		////! 三阶贝塞尔曲线的中间两个控制点的计算横跨四个节点，所以我们不需要担心B样条的中间节点是四重的
		////! 所以我们可以直接使用assert，但是首尾两个控制点横跨三个节点，他们有可能是一样的，
		////! 所以我们需要判断一下这三个节点是否一样，如果一样，说明这段曲线就是贝塞尔曲线，我们直接分解出来控制点即可
		////! 否则，这就是非均匀三次贝塞尔曲线，我们就需要计算中间控制点，可能会有重复计算的，但是这样可以统一了
		//
		//Vector2d p0_1 = p1;
		//assert(knotVec[j + 1] > knotVec[j - 2] && "Error: knot vector is not in ascending order.");
		//Vector2d p0_2 = ((knotVec[j] - knotVec[j - 2]) * controlPoints[j - 2] + (knotVec[j + 1] - knotVec[j]) * controlPoints[j - 3]) / (knotVec[j + 1] - knotVec[j - 2]);
		//assert(knotVec[j + 1] > knotVec[j - 1] && "Error: knot vector is not in ascending order.");
		//Vector2d p0 = ((knotVec[j] - knotVec[j - 1]) * p0_1 + (knotVec[j + 1] - knotVec[j]) * p0_2) / (knotVec[j + 1] - knotVec[j - 1]);

		//assert(knotVec[j + 3] > knotVec[j] && "Error: knot vector is not in ascending order.");
		//Vector2d p3_1 = ((knotVec[j + 1] - knotVec[j]) * controlPoints[j] + (knotVec[j + 3] - knotVec[j + 1]) * controlPoints[j - 1]) / (knotVec[j + 3] - knotVec[j]);
		//Vector2d p3_2 = p2;
		//assert(knotVec[j + 2] > knotVec[j] && "Error: knot vector is not in ascending order.");
		//Vector2d p3 = ((knotVec[j + 1] - knotVec[j]) * p3_1 + (knotVec[j + 2] - knotVec[j + 1]) * p3_2) / (knotVec[j + 2] - knotVec[j]);

		Eigen::Vector2d p0 = DeRecursion(knotVec, controlPoints, j - 1, 2, j);
		Eigen::Vector2d p1 = DeRecursion(knotVec, controlPoints, j - 1, 1, j);
		Eigen::Vector2d p2 = DeRecursion(knotVec, controlPoints, j - 1, 1, j + 1);
		Eigen::Vector2d p3 = DeRecursion(knotVec, controlPoints, j, 2, j + 1);

		std::vector<Vector2d> control_points = { p0, p1, p2, p3 };
		CubicBezier bezier(control_points, range);
		bezierVec.push_back(bezier);
	}
}

Eigen::Vector2d CubicBezier::DeRecursion(const vector<double>& knotVec, const vector<Vector2d>& controlPoints, const int& _lower, const int& _upper, const int& index)
{
	Vector2d ans;
	//int index = std::distance(knotVec.begin(), std::lower_bound(knotVec.begin(), knotVec.end(), t));
	if (!_upper && index - _lower <= 3 && index - _lower >= 0)
	{
		ans = controlPoints[_lower];
	}
	else if(_lower<=3 && index - _lower <= 3-_upper && index - _lower >= 0)
	{
		assert(knotVec[_lower - _upper + 4] > knotVec[_lower] && "Error: knot vector is not in ascending order.");
		double denominantor = knotVec[_lower - _upper + 4] - knotVec[_lower];
		ans = (knotVec[index] - knotVec[_lower]) * DeRecursion(knotVec, controlPoints, _lower, _upper - 1, index);
		ans += (knotVec[_lower - _upper + 4] - knotVec[index]) * DeRecursion(knotVec, controlPoints, _lower - 1, _upper - 1, index);
		ans /= denominantor;
	}
	return ans;
}

/*
*@brief: 输出三次贝塞尔曲线的控制点和有效域
****/
void operator<<(std::ostream& os, const CubicBezier& cb)
{
	for (auto& cp : cb.outContralPoints())
	{
		os<<cp.x()<<"\t"<<cp.y()<<"\t";
	}
	os << "Range: [ "<<cb.m_vecRange.first<<", "<<cb.m_vecRange.second<<" ]"  << std::endl;
}
