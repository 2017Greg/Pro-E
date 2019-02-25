#include "stdafx.h"
#include "creoRelatedAlgorithms.h"
#include <cmath>
#define PI 3.141592653

ProError VisitSolidSurface(ProSurface pSrf, ProError pStatus, ProAppData pData)
{
	SolidSurface* sSldSrf = (SolidSurface*)pData;
	ProArrayObjectAdd((ProArray*)&sSldSrf->pSrf  // 曲面向量 
		, -1                                     // 从尾部插入
		, 1                                      // 插入一个平面
		, &pSrf);                                // 平面的来源
	return PRO_TK_NO_ERROR;
}
ProError GetSolidSurface(ProMdl pSldMdl, std::vector<struct myplane> &allPlanes)
{
	ProSrftype pSrfType = PRO_SRF_PLANE;
	ProError err;
	SolidSurface sSldSrf;
	err = ProArrayAlloc(0, sizeof(ProSurface), 1, (ProArray*)&sSldSrf.pSrf);
	if (err != PRO_TK_NO_ERROR)
		return err;
	err = ProSolidSurfaceVisit((ProSolid)pSldMdl   // 实体句柄 
		, (ProSurfaceVisitAction)VisitSolidSurface     //  访问函数  如果它不返回 PRO_TK_NO_ERROR， 那么访问停止
		, NULL                   // 筛选函数 如果为空，表示访问所有的曲面 
		, &sSldSrf);             // 传递给访问函数和筛选函数的变量，用于存储访问中入队的实体表面
	if (err != PRO_TK_NO_ERROR)
		return err;
	int nSrfNum = 0;
	// 遍历实体表面 
	err = ProArraySizeGet(sSldSrf.pSrf, &nSrfNum);
	if (err != PRO_TK_NO_ERROR)
		return err;	

	for (int i = 0; i < nSrfNum; i++)
	{
		ProSurface pSrf = sSldSrf.pSrf[i];

		if (pSrfType != PRO_SRF_NONE)
		{
			ProGeomitem pSrfGeom;
			err = ProSurfaceToGeomitem((ProSolid)pSldMdl, pSrf, &pSrfGeom);
			if (err != PRO_TK_NO_ERROR)
				continue;
			ProGeomitemdata* pSrfGeomData;
			err = ProGeomitemdataGet(&pSrfGeom, &pSrfGeomData);

			if (err != PRO_TK_NO_ERROR)
				continue;
			ProSrftype pSrfType2 = pSrfGeomData->data.p_surface_data->type;

			if (pSrfType2 != pSrfType)
				continue;
			else {
				// 写出平面法向量的三个分量
				double a = pSrfGeomData->data.p_surface_data->srf_shape.plane.e3[0];
				double b = pSrfGeomData->data.p_surface_data->srf_shape.plane.e3[1];
				double c = pSrfGeomData->data.p_surface_data->srf_shape.plane.e3[2];
				// 写出平面的平移分量
				double x0 = pSrfGeomData->data.p_surface_data->srf_shape.plane.origin[0];
				double y0 = pSrfGeomData->data.p_surface_data->srf_shape.plane.origin[1];
				double z0 = pSrfGeomData->data.p_surface_data->srf_shape.plane.origin[2];
				// 写出平面的截距
				double d = -a*x0 - b*y0 - c*z0;
				// 则平面的方程为：ax+by+cz+d=0;
				struct myplane p;
				p.a = a;
				p.b = b;
				p.c = c;
				p.d = d;
				//process_plane(p);
				allPlanes.push_back(p);
			}
			ProGeomitemdataFree(&pSrfGeomData);
		}
	}

	ProArrayFree((ProArray*)&sSldSrf.pSrf);
	return err;
}


ProError annotationSelect(double &dimNorValue, double &upper_limit, double &lower_limit, std::string &annotype) {
	ProError err;
	ProMdl pSldMdl;
	err = ProMdlCurrentGet(&pSldMdl);
	if (PRO_TK_NO_ERROR != err )
	{
		AfxMessageBox(_T("please open a model first!"));
		return err;
	}
	/* 选择尺寸 */
	ProSelection *dims;
	int nDim = 0;
	err = ProSelect("annot_elem,feature ", 1, NULL, NULL, NULL, NULL, &dims, &nDim);   // gtol为 Geometric tolerance 
	if (PRO_TK_NO_ERROR != err || 1 != nDim)
	{
		AfxMessageBox(_T("errors occured in ProSelect()"));
		return err;
	}
	ProAnnotationElem annotationElem;
	err = ProSelectionModelitemGet(dims[0], &annotationElem);
	if (PRO_TK_NO_ERROR != err || 1 != nDim)
	{
		AfxMessageBox(_T("errors occured in ProSelectionModelitemGet()"));
		return err;
	}
	// annotation type
	ProAnnotationType type;
	ProAnnotationelemTypeGet(&annotationElem, &type);
	ProAnnotation annotation;
	ProAnnotationelemAnnotationGet(&annotationElem, &annotation);
	// 获取annotation的具体值：标称尺寸和公差，几何公差
	switch (type) {
	case PRO_ANNOT_TYPE_DRVDIM:
	{

		/* 根据ProDimension获取公称尺寸值 */
		ProDimensionNomvalueGet(&annotation, &dimNorValue);
		/*  根据ProDimension获取指定尺寸对象的公差  */
		ProDimensionToleranceGet(&annotation, &upper_limit, &lower_limit);
		/* 获取尺寸类型 */
		ProDimensiontype dim_type;
		ProDimensionTypeGet(&annotation, &dim_type);
		if (dim_type == PRODIMTYPE_LINEAR) {
			annotype = "linear_dimension";
			break;
		}
		if (dim_type == PRODIMTYPE_RADIUS) {
			annotype = "radius_dimension";
			break;
		}
		if (dim_type == PRODIMTYPE_DIAMETER) {
			annotype = "diameter_dimension";
			break;
		}
		if (!(dim_type == PRODIMTYPE_LINEAR || dim_type == PRODIMTYPE_RADIUS || dim_type == PRODIMTYPE_DIAMETER)) {
			annotype = "none";
			AfxMessageBox(_T("本检测系统暂无法处理该类型的尺寸"));
			break;
		}
	}
	break;
	case PRO_ANNOT_TYPE_GTOL:
	{
		wchar_t* value;
		ProGtolValueStringGet(&annotation, &value);
		dimNorValue = _wtof(value);
		/* 获取公差类型 */
		ProGtolType gtol_type;
		ProGtolTypeGet(&annotation, &gtol_type);
		if (gtol_type == PROGTOLTYPE_FLATNESS) {
			annotype = "gtol_flatness";
			break;
		}
		if (gtol_type == PROGTOLTYPE_CYLINDRICAL) {
			annotype = "gtol_cylindrical";
			break;
		}
		if (!(gtol_type == PROGTOLTYPE_FLATNESS || gtol_type == PROGTOLTYPE_CYLINDRICAL)) {
			annotype = "none";
			AfxMessageBox(_T("本检测系统暂无法处理该类型的公差"));
			break;
		}

	}
	break;
	default: AfxMessageBox(_T("您选择的注释对象不符合要求！！！请重新选择"));
		break;
	}
	return PRO_TK_NO_ERROR;
}
// FeatureGeomitemVisit的访问函数
ProError FeatureVisitAction(ProGeomitem* geomitem, ProError status, ProAppData data) {
	ProSurface surface;
	ProSrftype surf_type;
	status = ProGeomitemToSurface(geomitem, &surface);
	status = ProSurfaceTypeGet(surface, &surf_type);
	if (PRO_SRF_PLANE == surf_type || PRO_SRF_CYL == surf_type) {
		std::vector<ProGeomitem>* pf = (std::vector<ProGeomitem>*)data;
		pf->push_back(*geomitem);
	}
	return PRO_TK_NO_ERROR;
}

double ratio(double x, double y) {
	if (y == 0.0&&x != y)
		return 999999999;
	else if (y == 0.0&&x == 0.0) {
		return 1;
	}
	else {
		return x / y;
	}
}

ProError surfaceEquation(std::vector<myplane> &allPlanesOfUpperFeat
	, std::vector<mycylinder> &allCylindersOfUpperFeat
	, myplane &plane1
	, mycylinder &cylinder1)
{
	ProError err;
	// 选择面
	ProSelection *sels;
	int nSel = 0;
	err = ProSelect("surface", 1, NULL, NULL, NULL, NULL, &sels, &nSel);
	if (PRO_TK_NO_ERROR != err || 1 != nSel)
	{
		return err;
	}
	// 获取选择的面
	ProGeomitem geomSurface;
	err = ProSelectionModelitemGet(sels[0], &geomSurface);
	ProSurface surface;
	err = ProGeomitemToSurface(&geomSurface, &surface);
	ProGeomitemdata* surf_data;
	err = ProSurfaceDataGet(surface, &surf_data);

	// 定义存储平面数据的变量  
	ProSrftype srfType;   // 保存曲面类型
	ProUvParam uv_min;    // 保存极小值
	ProUvParam uv_max;    // 保存极大值
	ProSurfaceshapedata surf_shape;     // 保存平面的类型和数据
	int surf_id;
	ProSurfacedataGet((ProSurfacedata*)((surf_data->data).p_surface_data),
		&srfType, 
		uv_min, 
		uv_max, 
		NULL, 
		&surf_shape, 
		&surf_id);

	// 根据面的类型获取具体的数据
	/*CStringW cstrInfo;*/
	switch (srfType) {
	case PRO_SRF_PLANE:
	{
		//origin
		double x = surf_shape.plane.origin[0];
		double y = surf_shape.plane.origin[1];
		double z = surf_shape.plane.origin[2];
		// e1
		double i1 = surf_shape.plane.e1[0];
		double j1 = surf_shape.plane.e1[1];
		double k1 = surf_shape.plane.e1[2];
		// e2
		double i2 = surf_shape.plane.e2[0];
		double j2 = surf_shape.plane.e2[1];
		double k2 = surf_shape.plane.e2[2];
		// e3
		double i3 = surf_shape.plane.e3[0];
		double j3 = surf_shape.plane.e3[1];
		double k3 = surf_shape.plane.e3[2];
		// uvmin,uvmax
		double umin = uv_min[0];
		double vmin = uv_min[1];
		double umax = uv_max[0];
		double vmax = uv_max[1];
		// 存储手选平面
		plane1.a = i3;
		plane1.b = j3;
		plane1.c = k3;
		plane1.d = -x*i3 - y*j3 - z*k3;
		//process_plane(plane1);
		// 存储垂直平面
		double a[4] = { i1,i1,i2,i2 };
		double b[4] = { j1,j1,j2,j2 };
		double c[4] = { k1,k1,k2,k2 };
		double d[4];

		ProVector points[4];
		ProUvParam uv_points[4];

		uv_points[0][0] = umin;
		uv_points[0][1] = vmin;
		uv_points[1][0] = umax;
		uv_points[1][1] = vmin;
		uv_points[2][0] = umin;
		uv_points[2][1] = vmin;
		uv_points[3][0] = umin;
		uv_points[3][1] = vmax;

		ProSurfaceXyzdataEval(surface, uv_points[0], points[0], NULL, NULL, NULL);
		ProSurfaceXyzdataEval(surface, uv_points[1], points[1], NULL, NULL, NULL);
		ProSurfaceXyzdataEval(surface, uv_points[2], points[2], NULL, NULL, NULL);
		ProSurfaceXyzdataEval(surface, uv_points[3], points[3], NULL, NULL, NULL);

		d[0] = -a[0] * (points[0][0]) - b[0] * (points[0][1]) - c[0] * (points[0][2]);
		d[1] = -a[1] * (points[1][0]) - b[1] * (points[1][1]) - c[1] * (points[1][2]);
		d[2] = -a[2] * (points[2][0]) - b[2] * (points[2][1]) - c[2] * (points[2][2]);
		d[3] = -a[3] * (points[3][0]) - b[3] * (points[3][1]) - c[3] * (points[3][2]);

		// 存储垂直包围的四个平面
		for (int i = 0; i < 4; i++) {
			myplane plane0;
			plane0.a = a[i];
			plane0.b = b[i];
			plane0.c = c[i];
			plane0.d = d[i];
			//process_plane(plane0);
			allPlanesOfUpperFeat.push_back(plane0);
		}
	}
	break;
	case PRO_SRF_CYL:
	{
		// 圆柱半径
		cylinder1.radius =surf_shape.cylinder.radius;
		// 圆柱原点
		cylinder1.origin[0] = surf_shape.cylinder.origin[0];
		cylinder1.origin[1] = surf_shape.cylinder.origin[1];
		cylinder1.origin[2] = surf_shape.cylinder.origin[2];
		// 圆柱法向量
		cylinder1.e[0] = surf_shape.cylinder.e3[0];
		cylinder1.e[1] = surf_shape.cylinder.e3[1];
		cylinder1.e[2] = surf_shape.cylinder.e3[2];
		// uvmin,uvmax
		double umin = uv_min[0];
		double vmin = uv_min[1];
		double umax = uv_max[0];
		double vmax = uv_max[1];
		// 上下底面  
		/* 通过归纳法得出如下经验结论
		* 上下底面中心点的求解方法如下
		*  1、底面1： point1(x,y,z) = vmin * e3(x,y,z) + origin(x,y,z);
		*  2、底面2： point2(x,y,z) = vmax * e3(x,y,z) + origin(x,y,z);
		*/
		ProVector points[2];
		points[0][0] = cylinder1.origin[0] + vmin * cylinder1.e[0];
		points[0][1] = cylinder1.origin[1] + vmin * cylinder1.e[1];
		points[0][2] = cylinder1.origin[2] + vmin * cylinder1.e[2];
		// 由一个底面推导另一个底面
		//double t = (vmax - vmin) / pow(cylinder1.e[0] * cylinder1.e[0] + cylinder1.e[1] * cylinder1.e[1] + cylinder1.e[2] * cylinder1.e[2], 0.5);
		points[1][0] = cylinder1.origin[0] + vmax * cylinder1.e[0];
		points[1][1] = cylinder1.origin[1] + vmax * cylinder1.e[1];
		points[1][2] = cylinder1.origin[2] + vmax * cylinder1.e[2]; 

		double a[2] = { cylinder1.e[0],cylinder1.e[0] };
		double b[2] = { cylinder1.e[1],cylinder1.e[1] };
		double c[2] = { cylinder1.e[2],cylinder1.e[2] };
		double d[2];
		d[0] = -a[0] * points[0][0] - b[0] * points[0][1] - c[0] * points[0][2];
		d[1] = -a[1] * points[1][0] - b[1] * points[1][1] - c[1] * points[1][2];
		for (int i = 0; i < 2; i++) {
			myplane plane0;
			plane0.a = a[i];
			plane0.b = b[i];
			plane0.c = c[i];
			plane0.d = d[i];
			//process_plane(plane0);
			allPlanesOfUpperFeat.push_back(plane0);
		}
	}
	break;
	default:
		AfxMessageBox(_T("选取的面的类型不是平面或者圆柱面，暂不处理"));
		break;
	}
	return PRO_TK_NO_ERROR;
}

ProError getModelViewNames(ProMdl model, std::vector<std::string>& str_view_names) {    // bug exist -wk
	ProError status;
	ProLine* view_names;
	int numberOfNames;
	status = ProViewNamesGet(model,&view_names,NULL,&numberOfNames);
	if (status != PRO_TK_NO_ERROR) {
		AfxMessageBox(_T("获取视图名失败！ "));
		return status;
	}
	// 格式转换
	for (int i = 0; i < numberOfNames;i++) {
		int pSize = WideCharToMultiByte(CP_ACP, 0, view_names[i], wcslen(view_names[i]), NULL, 0, NULL, NULL);
		char* pCStrKey = new char[pSize+1];
		WideCharToMultiByte(CP_ACP, 0, view_names[i], wcslen(view_names[i])+1, pCStrKey, pSize+1, NULL, NULL);
		
		std::string str_temp = pCStrKey;
		delete[] pCStrKey;
		pCStrKey = NULL;
		str_view_names.push_back(str_temp);    
	}
	return status;
}

// 显示最佳视图 
void showView(ProName viewName)
{
	//Cdialog dlg;
	//dlg.DoModal();
	int status;
	ProMdl model;
	ProMdlCurrentGet(&model);
	ProView myView;
	status = ProViewRetrieve(model, viewName, &myView);
	if (status != PRO_TK_NO_ERROR) {
		CString cstr;
		cstr.Format(_T("%d"), status);
		AfxMessageBox(_T("errors occured in ProViewRetrieve()   ") + cstr);
		return;
	}
	status = ProViewRefit(model, myView);
	if (status != PRO_TK_NO_ERROR) {
		CString cstr;
		cstr.Format(_T("%d"), status);
		AfxMessageBox(_T("errors occured in ProViewRefit()   ") + cstr);
		return;
	}
}
	double cylinderNess;    // 存储圆柱度
	int size = cloud->points.size();
	bool flag = FALSE;   // flag=0表示圆柱母线不与X,Y,Z中任何坐标轴平行
	if (isDoubleEquals0(cylinder.i) && isDoubleEquals0(cylinder.j) && isDoubleEquals0(cylinder.k)) {
		AfxMessageBox(_T("提取的圆柱方程有误，无法对圆柱度误差进行评定！"));
		return -1.0;
	}
	// 当圆柱的母线与Z轴平行时
	if (isDoubleEquals0(cylinder.i) && isDoubleEquals0(cylinder.j)) {
		std::vector<double> ErrorsZ;
		for (int i = 0; i < size; i++) {
			double error;
			error = pow(pow((cloud->points[i].x - cylinder.x), 2.0) +
				pow(cloud->points[i].y - cylinder.y, 2.0), 0.5) - cylinder.r;
			ErrorsZ.push_back(error);
		}
		double minErrorZ = 9999999.0;
		double maxErrorZ = -9999999.0;
		for (int i = 0; i < size; i++) {
			if (minErrorZ > ErrorsZ[i])
				minErrorZ = ErrorsZ[i];
			if (maxErrorZ < ErrorsZ[i])
				maxErrorZ = ErrorsZ[i];
		}
		cylinderNess = maxErrorZ - minErrorZ;
		flag = TRUE;
	}

	// 当圆柱的母线与Y轴平行时
	if (isDoubleEquals0(cylinder.i) && isDoubleEquals0(cylinder.k)) {
		std::vector<double> ErrorsY;
		for (int i = 0; i < size; i++) {
			double error;
			error = pow(pow((cloud->points[i].x - cylinder.x), 2.0) +
				pow(cloud->points[i].z - cylinder.z, 2.0), 0.5) - cylinder.r;
			ErrorsY.push_back(error);
		}
		double minErrorY = 9999999.0;
		double maxErrorY = -9999999.0;
		for (int i = 0; i < size; i++) {
			if (minErrorY > ErrorsY[i])
				minErrorY = ErrorsY[i];
			if (maxErrorY < ErrorsY[i])
				maxErrorY = ErrorsY[i];
		}
		cylinderNess = maxErrorY - minErrorY;
		flag = TRUE;
	}

	// 当圆柱的母线与X轴平行时
	if (isDoubleEquals0(cylinder.j) && isDoubleEquals0(cylinder.k)) {
		std::vector<double> ErrorsX;
		for (int i = 0; i < size; i++) {
			double error;
			error = pow(pow((cloud->points[i].y - cylinder.y), 2.0) +
				pow(cloud->points[i].z - cylinder.z, 2.0), 0.5) - cylinder.r;
			ErrorsX.push_back(error);
		}
		double minErrorX = 9999999.0;
		double maxErrorX = -9999999.0;
		for (int i = 0; i < size; i++) {
			if (minErrorX > ErrorsX[i])
				minErrorX = ErrorsX[i];
			if (maxErrorX < ErrorsX[i])
				maxErrorX = ErrorsX[i];
		}
		cylinderNess = maxErrorX - minErrorX;
		flag = TRUE;
	}

	//  当圆柱母线不与任何坐标轴平行时  
	if(!flag){
		std::vector<double> Errors;
		// 将圆柱母线变换到与Z轴平行
		coordinateChange(cloud, cylinder);
		// 执行圆柱母线与Z轴平行的代码
		for (int i = 0; i < size; i++) {
			double error;
			error = pow(pow((cloud->points[i].x - cylinder.x), 2.0) +
				pow(cloud->points[i].y - cylinder.y, 2.0), 0.5) - cylinder.r;
			Errors.push_back(error);
		}
		double minErrorZ = 9999999.0;
		double maxErrorZ = -9999999.0;
		for (int i = 0; i < size; i++) {
			if (minErrorZ > Errors[i])
				minErrorZ = Errors[i];
			if (maxErrorZ < Errors[i])
				maxErrorZ = Errors[i];
		}
		cylinderNess = maxErrorZ - minErrorZ;
	}

	return cylinderNess;
}


void getBaryCenterOfSolid(ProMdl pSldMdl, std::vector<float> &BaryCenter) {
	ProError status;
	Pro3dPnt r_outline_points[2];
	ProSolidOutlineGet((ProSolid)pSldMdl, r_outline_points);     // 求解基坐标方向的外形轮廓
	for (size_t i = 0; i < 3;i++) {
		BaryCenter.push_back((r_outline_points[0][i] + r_outline_points[1][i])/2.0);
	}
}


void getOutLineOfSolid2(ProMdl pSldMdl, std::vector<float>& edges) {
	ProError status;
	Pro3dPnt r_outline_points[2];
	ProSolidOutlineGet((ProSolid)pSldMdl, r_outline_points);     // 求解基坐标方向的外形轮廓
	// 求边长并排序
	edges.resize(3);
	for (int i = 0; i < 3; i++) {
		edges[i] = r_outline_points[1][i] - r_outline_points[0][i];
	}
	// 插入排个序，从大到小
	for (int i = 1; i < 3; i++) {
		for (int j = i; j > 0; j--) {
			if (edges[j] > edges[j - 1]) {
				float temp = edges[j];
				edges[j] = edges[j - 1];
				edges[j - 1] = temp;
			}
		}
	}
}




