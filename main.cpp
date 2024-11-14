#include <iostream>

#include "curveFit.h"
#include <random>
#include "BezierCurve.h"
#include "BezierToBiArc.h"

static std::vector<Eigen::Vector2d> readTXTData(string filePath)
{
	std::vector<Eigen::Vector2d> ans;
	std::ifstream file(filePath);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << filePath << std::endl;
		return ans;
	}
	double x, y;
	while (file >> x >> y) {
		ans.emplace_back(Eigen::Vector2d(x, y));
	}
	file.close();
	return ans;
}



static std::vector<Eigen::Vector2d> readPLTData()
{
	std::vector<Eigen::Vector2d> ans;
	size_t size = 100;
	std::vector<std::vector<std::string>> vecString(size, std::vector<std::string>());
	int index = 0;
	std::string pltFilePath = "C:/Users/Zhushengb/source/repos/CurveFit/alphabet.plt";
	std::ifstream pltFile(pltFilePath);

	bool bBeginRead = false;
	bool bEndRead = false;

	if (!pltFile.is_open()) {
		std::cerr << "Error opening file: " << pltFilePath << std::endl;
		return ans;
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
				std::string x = str.substr(2, index - 2);
				std::string y = str.substr(index + 1, str.size() - index - 1);
				ans.emplace_back(Eigen::Vector2d(std::stod(x), std::stod(y)));
				//std::cout << "X : " << x << "\t" << "Y : " << y << std::endl;
			}
			//std::cout << str << std::endl;
		}
		num++;
		//std::cout << std::endl;
		//std::cout << std::endl;
		//std::cout << std::endl;

	}
	return ans;
}

int test1()
{
	const double PI = 3.14159265358979323846;
	double seed = std::chrono::system_clock::now().time_since_epoch().count();

	std::mt19937 gen(seed);
	std::normal_distribution<double> dis(0, 0.2);
	int n = 500;
	double X_MIN = 20.0;
	double X_MAX = 50.0;
	double X_RANGE = X_MAX - X_MIN;


	double radius = 500;
	Eigen::Vector2d circle_center(0.0, 0.0);

	//Eigen::Rotation2D<double> anticlockwise(PI / 2);
	//Eigen::Rotation2D<double> clockwise(-PI / 2);
	// 
	//for (int i = 0; i < n; i++)
	//{
	//	double y = 100 * std::sin(0.025 * i)+dis(gen);
	//	Eigen::Vector2d point(i, y);
	//	points.emplace_back(point);
	//}



	//for (int i = 0; i < n; i++)
	//{
	//	double x = X_MIN + X_RANGE * i / n;
	//	double y = 3e-2*std::pow(x-40,3)+4e-1*std::pow(x-40,2)-x;
	//	Eigen::Vector2d point(x, y);
	//	points.emplace_back(point);
	//}


	//for (int i = 0; i < n; i++)
	//{
	//	double t = i * 1.0 / n;
	//	double x = 10*std::pow(1-t,4)+40*t*std::pow(1-t,3)+180*std::pow(t,2)*std::pow(1-t,2)+240*(1-t)*std::pow(t,3)+70*std::pow(t,4);
	//	double y = 80 * t * std::pow(1 - t, 3) + 300 * std::pow(t, 2) * std::pow(1 - t, 2) + 80 * (1 - t) * std::pow(t, 3);
	//	Eigen::Vector2d point(x, y);
	//	points.emplace_back(point);
	//}

	//Eigen::Vector2d normal_vec(0, 1);
	//Eigen::Rotation2D<double> start_rotation(start_angle);
	//Eigen::Rotation2D<double> end_rotation(end_angle);
	//Eigen::Vector2d start_normal_vec = start_rotation * normal_vec;
	//Eigen::Vector2d end_normal_vec = end_rotation * normal_vec;

	double start_angle = 0;
	double end_angle = PI;

	string str = "C:\\Users\\Zhushengb\\source\\repos\\Xu_BiArcFitting\\data\\";

	auto startTime = std::chrono::high_resolution_clock::now();
	for (int index = 1; index <= 17; index++)
	{
		string fileName = str + std::to_string(index) + ".txt";
		std::vector<Eigen::Vector2d> points;
		points = readTXTData(fileName);
		if (points.empty())
		{
			continue;
		}


		BezierCurve bezier_curve(points);
		std::vector<Eigen::Vector2d> bezier_points = bezier_curve.outContralPoints();
		std::ofstream openfile("C:\\Users\\Zhushengb\\source\\repos\\Xu_BiArcFitting\\data\\bezier_curve.txt", std::ios::out | std::ios::app);
		if (openfile.is_open())
		{
			int index = 0;
			for (size_t i = 0; i < bezier_points.size(); i++)
			{
				++index;
				openfile << bezier_points[i].x() << "\t" << bezier_points[i].y() << "\t";
				if (index % 4 == 0)
				{
					openfile << std::endl;
				}
			}

		}
		else
		{
			std::cout << "open file failed" << std::endl;
		}
		openfile.close();


	}



	//double x = circle_center.x() + radius * std::cos(start_angle);
	//double y = circle_center.y() + radius * std::sin(start_angle);
	//Eigen::Vector2d start_point(x, y);
	//points.emplace_back(start_point);

	//for (int i = 1; i <= n; i++)
	//{
	//	double angle = start_angle + i * (end_angle - start_angle) / n;
	//	std::cout << angle << std::endl;
	//	double x = circle_center.x() + radius * std::cos(angle);
	//	double y = circle_center.y() + radius * std::sin(angle);
	//	Eigen::Vector2d point(x, y);
	//	points.emplace_back(point);
	//}



	auto endTime = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
	std::cout << " Spend time: " << duration.count() << " milliseconds" << std::endl;


	return 0;
}


int test2()
{
	std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());
	std::normal_distribution<double> dis(0, 200.0);
	std::vector<Eigen::Vector2d> contral_points;
	static const int NUM_OF_BEZIER_CURVE = 2;
	//for (int i = 0; i < NUM_OF_BEZIER_CURVE; i++)
	//{
	//	for (int j = 0; j < 4; j++)
	//	{
	//		double x = dis(gen);
	//		double y = dis(gen);
	//		Eigen::Vector2d point(x, y);
	//		contral_points.emplace_back(point);
	//	}
	//}

	contral_points.emplace_back(Eigen::Vector2d(500, 0));
	contral_points.emplace_back(Eigen::Vector2d(800, 400));
	contral_points.emplace_back(Eigen::Vector2d(400, 800));
	contral_points.emplace_back(Eigen::Vector2d(0,500));


	//std::fstream openfile("C:\\Users\\Zhushengb\\source\\repos\\Xu_BiArcFitting\\data\\bezier_curve.txt", std::ios::out | std::ios::app);
	//if (openfile.is_open())
	//{
	//	int index = 0;
	//	for (size_t i = 0; i < contral_points.size(); i++)
	//	{
	//		++index;
	//		openfile << contral_points[i].x() << "\t" << contral_points[i].y() << "\t";
	//		if (index % 4 == 0)
	//		{
	//			openfile << std::endl;
	//		}
	//	}
	//}

	auto startTime = std::chrono::high_resolution_clock::now();

	BezierToBiArc bezier_to_biarc(contral_points);
	bezier_to_biarc.fromCurvesToBiArc();
	std::vector<BiArc> biarcs = bezier_to_biarc.outBiArcs();

	std::fstream file("C:\\Users\\Zhushengb\\source\\repos\\Xu_BiArcFitting\\data\\biarc_curve.txt", std::ios::out | std::ios::app);
	if (file.is_open())
	{
		for (const auto& biarc : biarcs)
		{
			file << biarc;
		}
		file.close();

	}
	
	auto endTime = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
	std::cout << " Spend time: " << duration.count() << " milliseconds" << std::endl;

	return 0;
}


int main()
{
	std::vector<Eigen::Vector2d> contralPoints;
	CurveFit curveFit;
	curveFit.getData();
	curveFit.denseData();
	curveFit.doLineDetect();
	curveFit.doBiArcFit();
	vector<vector<Eigen::Vector2d>> curves = curveFit.getCurves();
	auto startTime = std::chrono::high_resolution_clock::now();
	for (auto& curve : curves)
	{
		BezierCurve bezier_curve(curve);
		vector<Eigen::Vector2d> contralpoint = bezier_curve.outContralPoints();
		contralPoints.insert(contralPoints.end(), contralpoint.begin(), contralpoint.end());
	}
	auto endTime = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
	std::cout << " Spend time: " << duration.count() << " milliseconds" << std::endl;

	{
		std::ofstream openfile("C:\\Users\\Zhushengb\\source\\repos\\Xu_BiArcFitting\\data\\bezier_curve.txt", std::ios::out | std::ios::app);
		if (openfile.is_open())
		{
			int index = 0;
			for (size_t i = 0; i < contralPoints.size(); i++)
			{
				++index;
				openfile << contralPoints[i].x() << "\t" << contralPoints[i].y() << "\t";
				if (index % 4 == 0)
				{
					openfile << std::endl;
				}
			}

		}
		else
		{
			std::cout << "open file failed" << std::endl;
		}
		openfile.close();
	}

	return 0;
}