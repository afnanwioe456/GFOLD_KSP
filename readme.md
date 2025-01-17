# GFOLD_KSP

## 结构

demo3_gfold.py：连游戏跑

GFOLD_direct_exec.py：建立求解模型

GFOLD_run.py：封装了求解器solver类，利用自带的主函数可以调用算法模拟运行并画一个图

EvilPlotting.py：画图相关代码，被GFOLD_run调用，自己不能跑

可调参数在params.txt中

## 主要改动

1.代码现代化改造，python3.6 -> python3.12，代码规范更新，各种库的版本核对，确保代码可以work；

2.取消画图，因为Unity版本问题经常画不上，想试试看的话可以把params.txt最后的debug_lines改成True；

3.新版cvxpy已经自己会使用cpp后端先编译问题再自动求解了，就移除了原仓库中关于代码生成的部分。

## Todo

1.控制算法有待改进，现在状态大体可用，但火箭下降到一定高度时竟然会上升，导致采样点失效，后期准备将引导算法改为顺序获取目标，飞到一个点附近强行切下一个点；

2.将项目接入基于PyQt的GUI中，绘制出凸优化求得的路径，并实时显示火箭路径，以便观察火箭实际运动和期望路径的关系，方便调试。

## Requirements

python 3.12+

通信部分：
krpc==0.5.4（别忘了游戏里mod版本的对应，不对应的话protobuf解析字节的时候协议对不上）
protobuf==5.29.3

算法部分：
numpy==2.2.1
scipy==1.15.1
cvxpy==1.6.0
ecos==2.0.14

简单起见直接用requirements装也可以(pip install -r requirements.txt)

## Reference

1.如何调试/修改的介绍视频链接：[【KSP】凸优化/GFOLD算法实现 代码环境搭建与排错](https://www.bilibili.com/video/BV1gN411p7Dm)；

2.凸优化参考论文：[Lossless Convexification of Nonconvex Control Bound and Pointing Constraints of the Soft Landing Optimal Control Problem](http://www.larsblackmore.com/iee_tcst13.pdf)

## License

此项目由[GFOLD_KSP](https://github.com/xdedss/GFOLD_KSP)改进而成，遵循GPLv3
