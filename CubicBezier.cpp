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
			vector<CubicBezier>& temp = transSplineToBezier(infile);
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
vector<CubicBezier>& CubicBezier::transSplineToBezier(std::ifstream& infile)
{
	static const string ENDLOG = "SPLINE";
	int DEGREE = 0;
	int KNOTSNUM = 0;
	int CONTROLNUM = 0;
	//! 读取到AcDbSpline字符串以后就开始读取样条曲线的各类数据了
	static auto ignore_N_line = [&infile](const int n)->void
		{
			for (int i = 0; i < n && infile.good(); ++i)
			{
				infile.ignore(50, '\n');
			}
	};
	//! 读取组码71、72、73所表示的样条曲线的次数、节点向量的长度、控制点的数量
	static auto get_param = [&infile, &DEGREE, &KNOTSNUM, &CONTROLNUM]() -> void
		{
			infile.ignore(50, '\n');
			string line;
			std::getline(infile, line);
			DEGREE = std::stoi(line);
			infile.ignore(50, '\n');
			std::getline(infile, line);
			KNOTSNUM = std::stoi(line);
			infile.ignore(50, '\n');
			std::getline(infile, line);
			CONTROLNUM = std::stoi(line);
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
	ignore_N_line(6);

	//! 依据组码72、73的值读取节点向量、控制点坐标
	{
		string line;
		//! 读取节点向量
		for (int i = 0; i < knots.size(); ++i)
		{
			ignore_N_line(1);
			std::getline(infile, line);
			knots[i] = std::stod(line);
		}


		{
			std::ofstream outfile("C:\\Users\\Zhushengb\\Desktop\\vector_graphic\\knots_011.txt", std::ios::app|std::ios::out);
			if (outfile.is_open())
			{
				for (int i = 1; i < knots.size(); i++)
				{
					outfile <<"TIME DIFFERENCRE: " <<i<<"\t" << knots[i] << "\t" << knots[i - 1] <<"\t is \t" << knots[i] - knots[i - 1] << std::endl;
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
	//! 将读取到的节点向量、控制点坐标转换为三次贝塞尔曲线
	vector<CubicBezier> ans(KNOTSNUM, CubicBezier());

	//for (int i = 0; i < KNOTSNUM-1; ++i)
	//{
	//	vector<Vector2d> subControlPoints(controlPoints.begin() + 3 * i, controlPoints.begin() + 3 * i + 4);
	//	CubicBezier cubicBezier(subControlPoints, std::make_pair(knots[i], knots[i+1]));
	//	ans[i] = cubicBezier;
	//}

	return ans;
}

/*
*@brief: 输出三次贝塞尔曲线的控制点和有效域
****/
void operator<<(std::ostream& os, const CubicBezier& cb)
{
	os << "Control Points: " << std::endl;
	for (auto& cp : cb.outContralPoints())
	{
		os<<cp<<std::endl;
	}
	os << "Range: [ "<<cb.m_vecRange.first<<", "<<cb.m_vecRange.second<<" ]"  << std::endl;
}
