#pragma once
#include<iostream>
#include<iomanip>
using namespace std;
struct OrbitRootNum
{
	double a;
	double e;
	double i;
	double w;
	double Omiga;
	double M;
	OrbitRootNum()//结构体初始化
	{
		a = 0;
		e = 0;
		i = 0;
		w = 0;
		Omiga = 0;
		M = 0;
	}
};
class rootnum
{
public:
	rootnum();//构造函数
	~rootnum();//析构函数
	bool SetRootNum(OrbitRootNum myOrbitRootNum);//输入轨道根数
	bool GetOrbitRootNum(OrbitRootNum& myOrbitRootNum, double t);//得到t时刻的轨道根数 
	bool GetOrbitRootNum(OrbitRootNum& myOrbitRootNum);//得到人卫单位下1Ts到100Ts的计算结果
private:
	double t;
	double j2;//j2项扁率摄动
	double n; //地球平均自转角速度
	double t0;//初始时刻
	double p; //常数，中间变量
	double r; //向径绝对距离
	double f; //真近点角
	double u ;//地球引力常数
	double rou;//角度和弧度转化
	double tao0;//天文时间单位
	double T;	//周期
	OrbitRootNum ORN;//用来承载中间量和返回结果，用这个来运算
	OrbitRootNum detaORN;//用来承载每项摄动影响的根数变化量
	OrbitRootNum Long1();//一阶长期项
	OrbitRootNum Long2();///二阶长期项
	OrbitRootNum ShortPeriod1();//一阶短周期项
	OrbitRootNum LongPeriod1();//一阶长周期项
	OrbitRootNum Add(OrbitRootNum OR1, OrbitRootNum OR2);//轨道根数相加
	OrbitRootNum Subtract(OrbitRootNum OR1, OrbitRootNum OR2);//轨道根数相减(OR1-OR2)
};

