#include "rootnum.h"
rootnum::rootnum()
{
	rou = 180 / acos(-1);
	j2 = 108263e-8;
	u = 398600.5e+9;
	tao0 = 806.8116347;//秒
	t0 = 0;
	t = 0;
	p = 0;
	n = 0;
	f = 0;
	r = 0;
	ORN.e = 0.1;
	ORN.i = 45 / rou;
	ORN.Omiga = 45 / rou;
	ORN.w = 45 / rou;
	ORN.M = 0;
	T = 120;//分钟
	ORN.a = pow(u * T * 60 * T * 60 / (4 * acos(-1) * acos(-1)), 1.0/3);//开普勒第三定律
}

rootnum::~rootnum()
{
}

bool rootnum::SetRootNum(OrbitRootNum myOrbitRootNum)
{
	ORN = myOrbitRootNum;
	return true;
}

bool rootnum::GetOrbitRootNum(OrbitRootNum& myOrbitRootNum, double t)//t为天文单位制
{
	OrbitRootNum detaORN;
	detaORN = ShortPeriod1();//计算一阶短周期项
	ORN = Subtract(ORN, detaORN);
	detaORN = LongPeriod1();//计算二阶长周期项
	ORN = Subtract(ORN, detaORN);//初始t0时刻平根数
	return true;
	t = t * tao0;
	detaORN = Long1();//计算一阶短周期项
	ORN = Add(ORN, detaORN);
	detaORN = Long2();//计算二阶短周期项
	ORN = Add(ORN, detaORN);
	n = sqrt(u / (ORN.a * ORN.a * ORN.a));
	ORN.M = ORN.M + n * (t - t0);//t时刻的平根数
	detaORN = Long1();//计算t时刻一阶短周期项
	ORN = Add(ORN, detaORN);
	detaORN = Long2();//计算t时刻二阶短周期项
	ORN = Add(ORN, detaORN);//t时刻的瞬时轨道根数
	while (ORN.M > 2 * acos(-1))
	{
		ORN.M -= 2 * acos(-1);
	}
	cout << "Ts=" << t << ":" << endl;
	cout << "a=" << ORN.a << "   " << "i=" << ORN.i << "   " << "e=" << ORN.e << "   " << "w=" << ORN.w << "   " << "Omiga=" << ORN.Omiga << "   " << "M=" << ORN.M << endl;
	myOrbitRootNum = ORN;
}

bool rootnum::GetOrbitRootNum(OrbitRootNum& myOrbitRootNum)
{
	OrbitRootNum detaORN;
	detaORN = ShortPeriod1();//计算一阶短周期项
	ORN = Subtract(ORN, detaORN);
	detaORN = LongPeriod1();//计算二阶长周期项
	ORN = Subtract(ORN, detaORN);//初始t0时刻平根数
	for (int i = 0; i < 100; i++)
	{
		t = (i + 1) * tao0;
		detaORN = Long1();//计算一阶短周期项
		ORN = Add(ORN, detaORN);
		detaORN = Long2();//计算二阶短周期项
		ORN = Add(ORN, detaORN);
		n = sqrt(u / (ORN.a * ORN.a * ORN.a));
		ORN.M = ORN.M + n * (t - t0);//t时刻的平根数
		detaORN = Long1();//计算t时刻一阶短周期项
		ORN = Add(ORN, detaORN);
		detaORN = Long2();//计算t时刻二阶短周期项
		ORN = Add(ORN, detaORN);//t时刻的瞬时轨道根数
		while (ORN.M > 2 * acos(-1))
		{
			ORN.M -= 2 * acos(-1);
		}
		cout << "Ts=" << i + 1 <<":"<< endl;
		//控制小数的输出位数
		cout << setiosflags(ios::fixed) << setprecision(15) << "a=" << ORN.a << "   " << "i=" << ORN.i << "   " << "e=" << ORN.e << "   " << "w=" << ORN.w << "   " << "Omiga=" << ORN.Omiga << "   " << "M=" << ORN.M << endl << endl;
		myOrbitRootNum = ORN;
	}
	return true;
}

OrbitRootNum rootnum::Long1()
{
	double E0, E;
	p = ORN.a * (1 - ORN.e * ORN.e);
	n = sqrt(u / (ORN.a * ORN.a * ORN.a));
	E0 = ORN.M;
	E = ORN.M + ORN.e * sin(E0);
	while (abs(E - E0) > 1e-8)
	{
		E0 = E;
		E = ORN.M + ORN.e * sin(E0);

	}
	f = atan2(sqrt(1 - ORN.e * ORN.e) * sin(E), cos(E) - ORN.e);
	if (f < 0)
		f = f + 2 * acos(-1);
	r = ORN.a * (1 - ORN.e * cos(E));
    OrbitRootNum detaORN;//用来承载每项摄动影响的根数变化量
	detaORN.a = 0;
	detaORN.e = 0;
	detaORN.i = 0;
	detaORN.Omiga = -3 * j2 * cos(ORN.i) * n * (t - t0) / (2 * p * p);
	detaORN.w = 3 * j2 * (2 - 5 * sin(ORN.i) * sin(ORN.i) / 2) * n / (2 * p * p) * (t - t0);
	detaORN.M = 3 * j2 * (1 - 3 * sin(ORN.i) * sin(ORN.i) / 2) * sqrt(1 - ORN.e * ORN.e) * n / (2 * p * p) * (t - t0);
    return detaORN;
}

OrbitRootNum rootnum::Long2()
{
	p = ORN.a * (1 - ORN.e * ORN.e);
	n = sqrt(u / (ORN.a * ORN.a * ORN.a));
	double E0, E;
	E0 = ORN.M;
	E = ORN.M + ORN.e * sin(E0);
	while (abs(E - E0) > 1e-8)
	{
		E0 = E;
		E = ORN.M + ORN.e * sin(E0);

	}
	f = atan2(sqrt(1 - ORN.e * ORN.e) * sin(E), cos(E) - ORN.e);
	if (f < 0)
		f = f + 2 * acos(-1);
	r = ORN.a * (1 - ORN.e * cos(E));
	OrbitRootNum detaORN;//用来承载每项摄动影响的根数变化量
	detaORN.a = 0;
	detaORN.e = 0;
	detaORN.i = 0;
	detaORN.Omiga = -pow(3 * j2 / (2 * p * p), 2) * cos(ORN.i) * ((3 / 2 + 1 * ORN.e * ORN.e / 6 + sqrt(1 - ORN.e * ORN.e)) - sin(ORN.i) * sin(ORN.i) * (5 / 3 - 5 / 24 * ORN.e * ORN.e + 3 / 2 * sqrt(1 - ORN.e * ORN.e))) * n * (t - t0);
	detaORN.w = pow(3 * j2 / (2 * p * p), 2) * (4 + 7 / 12 * ORN.e * ORN.e + 2 * sqrt(1 - ORN.e * ORN.e) - sin(ORN.i) * sin(ORN.i) * (130 / 12 + 3 * ORN.e * ORN.e / 8 + 11 * sqrt(1 - ORN.e * ORN.e) / 2) + pow(sin(ORN.i), 4) * (215 / 48 - 15 * ORN.e * ORN.e / 32 + 15 / 4 * sqrt(1 - ORN.e * ORN.e))) * n * (t - t0);
	detaORN.M = pow(3 * j2 / (2 * p * p), 2) * sqrt(1 - ORN.e * ORN.e) * (0.5 * pow(1 - 3 / 2 * sin(ORN.i) * sin(ORN.i), 2) * sqrt(1 - ORN.e * ORN.e) + 2.5 + 10 / 3 * ORN.e * ORN.e - sin(ORN.i) * sin(ORN.i) * (19 / 3 + 26 * ORN.e * ORN.e / 3) + pow(sin(ORN.i), 4) * (233 / 48 + 103 / 12 * ORN.e * ORN.e) + pow(ORN.e, 4) / (1 - ORN.e * ORN.e) * (35 / 12 - 35 * sin(ORN.i) * sin(ORN.i) / 4 + 315 * pow(sin(ORN.i), 4) / 32)) * n * (t - t0);
	return detaORN;
}

OrbitRootNum rootnum::ShortPeriod1()
{
	double E0, E;
	E0 = ORN.M;
	E = ORN.M + ORN.e * sin(E0);
	while (abs(E - E0) > 1e-8)
	{
		E0 = E;
		E = ORN.M + ORN.e * sin(E0);

	}
	f = atan2(sqrt(1 - ORN.e * ORN.e) * sin(E), cos(E) - ORN.e);
	if (f < 0)
		f = f + 2 * acos(-1);
	r = ORN.a * (1 - ORN.e * cos(E));
	p = ORN.a * (1 - ORN.e * ORN.e);
	n = sqrt(u / (ORN.a * ORN.a * ORN.a));
	OrbitRootNum detaORN;//用来承载每项摄动影响的根数变化量
	detaORN.a = j2 / ORN.a * ((1 - 3 * sin(ORN.i) * sin(ORN.i) / 2) * ((pow(ORN.a / r, 3)) - pow(1 - ORN.e * ORN.e, -2 / 3)) + 3 * sin(ORN.i) * sin(ORN.i) * pow(ORN.a / r, 3) * cos(2 * (ORN.w + f)));
	detaORN.e = j2 / (2 * ORN.a * ORN.a) * ((1 - ORN.e * ORN.e) / ORN.e) * ((1 - 3 * sin(ORN.i) * sin(ORN.i) / 2) * ((pow(ORN.a / r, 3)) - pow(1 - ORN.e * ORN.e, -2 / 3)) + 3 * sin(ORN.i) * sin(ORN.i) * pow(ORN.a / r, 3) * cos(2 * (ORN.w + f)) - 3 * sin(ORN.i) * sin(ORN.i) / (2 * pow(1 - ORN.e * ORN.e, 2)) * (ORN.e * cos(f + 2 * ORN.w)) + cos(2 * f + 2 * ORN.w) + ORN.e / 3 * cos(3 * f + 2 * ORN.w));
	detaORN.i = 3 * j2 * sin(2 * ORN.i) * (ORN.e * cos(f + 2 * ORN.w) * cos(2 * f + 2 * ORN.w) + ORN.e / 3 * cos(3 * f + 2 * ORN.w));
	detaORN.w = 3 * j2 / (2 * p * p) * ((2 - 5 * sin(ORN.i) * sin(ORN.i) / 2) * (f - ORN.M + ORN.e * sin(f)) + (1 - 3 * sin(ORN.i) * sin(ORN.i) / 2) * ((1 / ORN.e - ORN.e / 4) * sin(f) + sin(2 * f) / 2 + ORN.e / 12 * sin(3 * f)) - sin(f + 2 * ORN.w) * (sin(ORN.i) * sin(ORN.i) / 4 / ORN.e + (1 / 2 - 15 * sin(ORN.i) * sin(ORN.i) / 16) * ORN.e) - (1 / 2 - 5 * sin(ORN.i) * sin(ORN.i) / 4) * sin(2 * f + 2 * ORN.w) + (7 / 12 / ORN.e * sin(ORN.i) * sin(ORN.i) - (1 / 6 - 19 / 48 * sin(ORN.i) * sin(ORN.i)) * ORN.e) * sin(3 * f + 2 * ORN.w) + (3 / 8 * sin(ORN.i) * sin(ORN.i)) * sin(4 * f + 2 * ORN.w) + ORN.e / 16 * sin(ORN.i) * sin(ORN.i) * (sin(5 * f + 2 * ORN.w) + sin(f - 2 * ORN.w)));
	detaORN.Omiga = -3 * j2 * cos(ORN.i) / (2 * p * p) * ((f - ORN.M + ORN.e * sin(f)) - 0.5 * (ORN.e * sin(f + 2 * ORN.w) + sin(2 * f + 2 * ORN.w) + ORN.e * sin(3 * f + 2 * ORN.w) / 3));
	detaORN.M = sqrt(1 - ORN.e * ORN.e) * 3 * j2 / (2 * p * p) * (-(1 - 3 / 2 * sin(ORN.i) * sin(ORN.i)) * ((1 / ORN.e - ORN.e / 4) * sin(f) + sin(2 * f) / 2 + ORN.e * sin(3 * f) / 12) + ((1 / (4 * ORN.e) + 15 / 16 * ORN.e) * sin(ORN.i) * sin(ORN.i)) * sin(f + 2 * ORN.w) - ((7 / 12 / ORN.e - ORN.e / 48) * sin(ORN.i) * sin(ORN.i) * sin(3 * f + 2 * ORN.w) + 3 / 8 * sin(ORN.i) * sin(ORN.i) * sin(4 * f + 2 * ORN.w)) - ORN.e / 16 * sin(ORN.i) * sin(ORN.i) * (sin(5 * f + 2 * ORN.w) + sin(f - 2 * ORN.w)));
	return detaORN;
}

OrbitRootNum rootnum::LongPeriod1()
{
	double E0, E;
	E0 = ORN.M;
	E = ORN.M + ORN.e * sin(E0);
	while (abs(E - E0) > 1e-8)
	{
		E0 = E;
		E = ORN.M + ORN.e * sin(E0);

	}
	f = atan2(sqrt(1 - ORN.e * ORN.e) * sin(E), cos(E) - ORN.e);
	if (f < 0)
		f = f + 2 * acos(-1);
	p = ORN.a * (1 - ORN.e * ORN.e);
	n = sqrt(u / (ORN.a * ORN.a * ORN.a));
	r = ORN.a * (1 - ORN.e * cos(E));
	OrbitRootNum detaORN;//用来承载每项摄动影响的根数变化量
	//一阶长周期项
	detaORN.a = 0;
	detaORN.i = -(3 * j2) / (2 * p * p) * sin(2 * ORN.i) / (4 - 5 * sin(ORN.i) * sin(ORN.i)) * (7 / 24 - 5 * sin(ORN.i) * sin(ORN.i) / 16) * pow(ORN.e, 2) * cos(2 * ORN.w);
	detaORN.Omiga = -(3 * j2) / (2 * p * p) * cos(ORN.i) / pow(4 - 5 * sin(ORN.i) * sin(ORN.i), 2) * (7 / 3 - 5 * sin(ORN.i) * sin(ORN.i) + 25 / 8 * pow(sin(ORN.i), 4)) * pow(ORN.e, 2) * sin(2 * ORN.w);
	detaORN.e = -((1 - ORN.e * ORN.e) / pow(ORN.e, 2) * tan(ORN.i)) * detaORN.i;//可能有问题
	detaORN.w = -(3 * j2) / (2 * p * p) / pow(4 - 5 * sin(ORN.i) * sin(ORN.i), 2) * (pow(sin(ORN.i), 2) * (25 / 3 - 245 / 12 * sin(ORN.i) * sin(ORN.i) + 25 / 2 * pow(sin(ORN.i), 4)) - ORN.e * ORN.e * (7 / 3 - 17 / 2 * sin(ORN.i) * sin(ORN.i) + 65 / 6 * pow(sin(ORN.i), 4) - 75 / 16 * pow(sin(ORN.i), 6))) * sin(2 * ORN.w);
	detaORN.M = (3 * j2) / (2 * p * p) * sin(ORN.i) * sin(ORN.i) / pow(4 - 5 * sin(ORN.i) * sin(ORN.i), 2) * sqrt(1 - ORN.e * ORN.e) * ((25 / 3 - 245 / 12 * sin(ORN.i) * sin(ORN.i) + 25 / 2 * pow(sin(ORN.i), 4)) - ORN.e * ORN.e * (4 - 5 * sin(ORN.i) * sin(ORN.i)) * (7 / 12 - 5 / 8 * sin(ORN.i) * sin(ORN.i))) * sin(2 * ORN.w);
	return detaORN;
}

OrbitRootNum rootnum::Add(OrbitRootNum OR1, OrbitRootNum OR2)
{
	OrbitRootNum AddORN;
	AddORN.a = OR1.a + OR2.a;
	AddORN.e = OR1.e + OR2.e;
	AddORN.i = OR1.i + OR2.i;
	AddORN.w = OR1.w + OR2.w;
	AddORN.Omiga = OR1.Omiga + OR2.Omiga;
	AddORN.M = OR1.M + OR2.M;
	return AddORN;
}
OrbitRootNum rootnum::Subtract(OrbitRootNum OR1, OrbitRootNum OR2)
{
	OrbitRootNum AddORN;
	AddORN.a = OR1.a - OR2.a;
	AddORN.e = OR1.e - OR2.e;
	AddORN.i = OR1.i - OR2.i;
	AddORN.w = OR1.w - OR2.w;
	AddORN.Omiga = OR1.Omiga - OR2.Omiga;
	AddORN.M = OR1.M - OR2.M;
	return AddORN;
}