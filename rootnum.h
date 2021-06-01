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
	OrbitRootNum()//�ṹ���ʼ��
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
	rootnum();//���캯��
	~rootnum();//��������
	bool SetRootNum(OrbitRootNum myOrbitRootNum);//����������
	bool GetOrbitRootNum(OrbitRootNum& myOrbitRootNum, double t);//�õ�tʱ�̵Ĺ������ 
	bool GetOrbitRootNum(OrbitRootNum& myOrbitRootNum);//�õ�������λ��1Ts��100Ts�ļ�����
private:
	double t;
	double j2;//j2������㶯
	double n; //����ƽ����ת���ٶ�
	double t0;//��ʼʱ��
	double p; //�������м����
	double r; //�򾶾��Ծ���
	double f; //������
	double u ;//������������
	double rou;//�ǶȺͻ���ת��
	double tao0;//����ʱ�䵥λ
	double T;	//����
	OrbitRootNum ORN;//���������м����ͷ��ؽ���������������
	OrbitRootNum detaORN;//��������ÿ���㶯Ӱ��ĸ����仯��
	OrbitRootNum Long1();//һ�׳�����
	OrbitRootNum Long2();///���׳�����
	OrbitRootNum ShortPeriod1();//һ�׶�������
	OrbitRootNum LongPeriod1();//һ�׳�������
	OrbitRootNum Add(OrbitRootNum OR1, OrbitRootNum OR2);//����������
	OrbitRootNum Subtract(OrbitRootNum OR1, OrbitRootNum OR2);//����������(OR1-OR2)
};

