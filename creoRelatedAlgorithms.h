#pragma once
#include "stdafx.h"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>

// creo相关头文件 
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
	double d;          // 平面一般方程的四个变量（a,b,c）为方向向量，d为截距 
};

struct mycylinder {
	double radius;        // 圆柱半径
	ProPoint3d  origin;   // 原点
	ProPoint3d  e;        // 圆柱法向量
};

struct SolidSurface
{
	ProSurface* pSrf;
};
ProError VisitSolidSurface(ProSurface pSrf, ProError pStatus, ProAppData pData);
/*GetSolidSurface：提取当前模型的所有平面的所有方程，用于配准 
* 输入参数：ProMdl pSldMdl，CAD模型
* 输出参数：std::vector<struct myplane> &allPlanes， 从模型中提取的所有平面
*/
ProError GetSolidSurface(ProMdl pSldMdl, std::vector<struct myplane> &allPlanes);

/*annotationSelect：提取尺寸的标称值和上下公差/ 提取几何公差的值 
* 输出参数1：double &dimNorValue，存储尺寸的标称值 
* 输出参数2：double &upper_limit，存储尺寸的上公差
* 输出参数3：double &lower_limit，存储尺寸的下公差
* 输出参数4：std::string &annotype)，注释类型
*/
ProError annotationSelect(double &dimNorValue          
	, double &upper_limit                        
	, double &lower_limit                          
	, std::string &annotype); 

	
ProError FeatureVisitAction(ProGeomitem* geomitem, ProError status, ProAppData data);
/*surfaceEquation：提取尺寸关联几何以及尺寸上层特征的包络曲面，包括平面和圆柱面  
* 输出参数1： std::vector<myplane> &allPlanesOfUpperFeat，上层特征的所有平面
* 输出参数2： std::vector<mycylinder> &allCylinderOfupperFeat，上层特征的所有柱面
* 输出参数3： myplane &plane1，尺寸关联的平面
* 输出参数4： mycylinder &cylinder1)，尺寸关联的曲面
*/
ProError surfaceEquation(std::vector<myplane> &allPlanesOfUpperFeat   
	, std::vector<mycylinder> &allCylinderOfupperFeat           
	, myplane &plane1                                           
	, mycylinder &cylinder1);   
	                     
						 
/* getModelViewNames：获取模型的所有视图名 
* 输入参数：ProMdl model，CAD模型
* 输出参数：std::vector<std::string>& str_view_names， 模型中保存的所有视图
*/
ProError getModelViewNames(ProMdl model,std::vector<std::string>& str_view_names);

/* showView: 根据视图名显示视图
*输入参数：ProName viewName，视图名
*/
void showView(ProName viewName);

/* getBaryCenterOfSolid功能：获取零件的外包络长方体中心点
* 输入参数：ProMdl pSldMdl; proe格式模型
* 输出参数：std::vector<float> &BaryCenter; 重心坐标
*/
void getBaryCenterOfSolid(ProMdl pSldMdl, std::vector<float> &BaryCenter);

/* getOutLineOfSolid功能：获取零件的外包络长方体长宽高
* 输入参数：ProMdl pSldMdl; proe格式模型
* 输出参数：std::vector<float>& edges; 长宽高，按长度降序排列
*/
void getOutLineOfSolid2(ProMdl pSldMdl, std::vector<float>& edges);

