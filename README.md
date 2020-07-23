# Pro/E

## 描述/description  
Pro/E二次开发的部分代码


## 代码使用方法/code usage  

1 参考 https://blog.csdn.net/weixin_38354109/article/details/79397218 配置protoolkit开发环境  
2 将两个文件.h和.cpp文件添加到MFC-dll工程中  
3 根据功能要求调用相关的函数  
4 编译工程，生成dll文件，编写.dat配置文件和.txt菜单文件，用creo/proe加载工程，即可看到效果  


## 代码功能介绍
1 GetSolidSurface()  # 获取模型（model）的所有平面的方程，存储在 allplanes变量中  
2 annotationSelect() # 获取手动选择尺寸的标称值、上下公差值、注释类型  
3 surfaceEquation()  # 获取手选曲面(surface)的方程  
4 getModelViewNames() # 获取当前模型的所有视图名  
5 getBaryCenterOfSolid() # 获取CAD模型的几何中心  
6 getOutLineOfSolid2()   # 获取CAD模型的外接矩的长宽高  


