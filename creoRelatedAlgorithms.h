#pragma once
#include "stdafx.h"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>

// creo���ͷ�ļ� 
#include <ProSolid.h>
#include <ProSelection.h>
#include <ProDimension.h>
#include <ProAnnotation.h>
#include <ProUtil.h>
#include <ProMessage.h>
#include <ProWindows.h>
#include <ProMessage.h>


struct myplane {
	double a;
	double b;
	double c;
	double d;          // ƽ��һ�㷽�̵��ĸ�������a,b,c��Ϊ����������dΪ�ؾ� 
};

struct mycylinder {
	double radius;        // Բ���뾶
	ProPoint3d  origin;   // ԭ��
	ProPoint3d  e;        // Բ��������
};

struct SolidSurface
{
	ProSurface* pSrf;
};
ProError VisitSolidSurface(ProSurface pSrf, ProError pStatus, ProAppData pData);
/*GetSolidSurface����ȡ��ǰģ�͵�����ƽ������з��̣�������׼ 
* ���������ProMdl pSldMdl��CADģ��
* ���������std::vector<struct myplane> &allPlanes�� ��ģ������ȡ������ƽ��
*/
ProError GetSolidSurface(ProMdl pSldMdl, std::vector<struct myplane> &allPlanes);

/*annotationSelect����ȡ�ߴ�ı��ֵ�����¹���/ ��ȡ���ι����ֵ 
* �������1��double &dimNorValue���洢�ߴ�ı��ֵ 
* �������2��double &upper_limit���洢�ߴ���Ϲ���
* �������3��double &lower_limit���洢�ߴ���¹���
* �������4��std::string &annotype)��ע������
*/
ProError annotationSelect(double &dimNorValue          
	, double &upper_limit                        
	, double &lower_limit                          
	, std::string &annotype); 

	
ProError FeatureVisitAction(ProGeomitem* geomitem, ProError status, ProAppData data);
/*surfaceEquation����ȡ�ߴ���������Լ��ߴ��ϲ������İ������棬����ƽ���Բ����  
* �������1�� std::vector<myplane> &allPlanesOfUpperFeat���ϲ�����������ƽ��
* �������2�� std::vector<mycylinder> &allCylinderOfupperFeat���ϲ���������������
* �������3�� myplane &plane1���ߴ������ƽ��
* �������4�� mycylinder &cylinder1)���ߴ����������
*/
ProError surfaceEquation(std::vector<myplane> &allPlanesOfUpperFeat   
	, std::vector<mycylinder> &allCylinderOfupperFeat           
	, myplane &plane1                                           
	, mycylinder &cylinder1);   
	                     
						 
/* getModelViewNames����ȡģ�͵�������ͼ�� 
* ���������ProMdl model��CADģ��
* ���������std::vector<std::string>& str_view_names�� ģ���б����������ͼ
*/
ProError getModelViewNames(ProMdl model,std::vector<std::string>& str_view_names);

/* showView: ������ͼ����ʾ��ͼ
*���������ProName viewName����ͼ��
*/
void showView(ProName viewName);

/* getBaryCenterOfSolid���ܣ���ȡ���������糤�������ĵ�
* ���������ProMdl pSldMdl; proe��ʽģ��
* ���������std::vector<float> &BaryCenter; ��������
*/
void getBaryCenterOfSolid(ProMdl pSldMdl, std::vector<float> &BaryCenter);

/* getOutLineOfSolid���ܣ���ȡ���������糤���峤���
* ���������ProMdl pSldMdl; proe��ʽģ��
* ���������std::vector<float>& edges; ����ߣ������Ƚ�������
*/
void getOutLineOfSolid2(ProMdl pSldMdl, std::vector<float>& edges);

