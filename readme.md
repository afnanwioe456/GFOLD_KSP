# GFOLD_KSP

了解原仓库的情况请移步下一小节。

Mr.G进行的修改：

1.代码现代化改造，python3.6 -> python3.12，代码规范更新，各种库的版本核对，确保代码可以work；

2.暂时仅确保直接执行work，代码生成的部分后面再说；

3.取消画图，因为Unity版本问题经常画不上，想试试看的话可以把params.txt最后的debug_lines改成True；

4.如何调试/修改的介绍视频链接：https://www.bilibili.com/video/BV1gN411p7Dm ；

5.

# 以下内容来自原仓库的readme

连接krpc并调用求解器进行计算和自动降落的主要程序是demo3_gfold.py

GFOLD_codegen.py - 生成求解的c代码

GFOLD_direct_exec.py - 直接cvxpy求解

GFOLD_run.py - 封装了solver类，可以调用上述两个方法求解

可调参数在params.txt中

## Requirements

python 3.x （推荐3.6）

setuptools 版本**小于等于** 57.5.0（否则后续安装cvxpy可能出现错误error in cvxpy setup command: use_2to3 is invalid.）

scipy==1.2.1，CVXcanon==0.1.0（这两个是cvxpy的依赖项，其当前最新版已不兼容我们使用的cvxpy版本，所以需要手动安装）

cvxpy==0.4.11（为了兼容cvxpy_codegen，不能使用新的1.x版本）

[cvxpy_codegen](https://github.com/moehle/cvxpy_codegen) （注：在python3.x中可能需要把setup.py中"python >= 2.7"一行删去才能正常安装）

[krpc](http://krpc.github.io/krpc/)

## Codegen

运行GFOLD_codegen.py可以生成效率较高的c代码并作为python包安装到python环境里。

注意生成的代码可能有许多错误需要看情况修正，因此我写了postprocess.py来自动改正一些常见问题，但是不保证在所有环境下都能完美工作。

上述整个代码生成+改错+编译安装过程可以通过运行build.bat自动完成，如果一切正常，运行完之后python里就会多出两个包 gfold_solver_p3 和 gfold_solver_p4 

## 直接在python中求解

如果要直接cvxpy求解（耗时可能达到c的10倍以上）请把params.txt中direct改成True，然后运行demo3_gfold.py

## 测试

运行python ./GFOLD_run.py，可以用默认的测试数据+codegen的代码求解并显示结果（必须要先完成codegen的步骤）

加参数运行python ./GFOLD_run.py direct，则是同样的默认数据用直接求解的方法求解并显示结果

如果默认数据能够得出正确结果（不会infeasible），说明求解器没有问题。

## License

此项目由[G-FOLD-Python](https://github.com/jonnyhyman/G-FOLD-Python)改进而成，遵循GPLv3
