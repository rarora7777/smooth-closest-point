------ Compilation

Before computing the Euclidean projection using the MATLAB code, first compile decimator:
cd decimator
cl decimator.cc /I OPENMESH_ROOT/include /link /LIBPATH:OPENMESH_ROOT/lib OpenMeshCore.lib OpenMeshTools.lib

Compile Phong projection code to a dynamic library (DLL) using
cd src/phong
cl /LD Phong.cpp Trianglephong.cpp /I "../utils" /I ".." /I EIGEN3_ROOT /I LIBIGL_ROOT/include

------ Introduction

Reference implementation for the paper Weighted Averages on Surfaces
Daniele Panozzo, Ilya Baran, Olga Diamanti, Olga Sorkine
SIGGRAPH 2013

This code is free to use for non-commercial applications, please
cite our paper if you use it in your project.

@article{Panozzo:WA:2013,
author = {Daniele Panozzo and Ilya Baran and Olga Diamanti and Olga Sorkine-Hornung},
title = {Weighted Averages on Surfaces},
journal = {ACM Transactions on Graphics (proceedings of ACM SIGGRAPH)},
volume = {32},
number = {4},
year = {2013},
}

------ Content of the package:

The code is divided into three parts:

1) Euclidean embedding metric computation (matlab) Section 3.5

This can be used independently, the script that contains the implementation is WA_precompute.m.
It depends on the fast_marching toolbox from Gabriel Peyre (http://www.mathworks.com/matlabcentral/fileexchange/6110-toolbox-fast-marching) and on OpenMesh(http://www.openmesh.org).

2) Phong Projection (C++)

This can be used independently, it is contained two files (TrianglePhong.h and TrianglePhong.cpp) and it only requires Eigen to be compiled. See the comments in the header file for more details.

3) Weighted Averages (C++)

It depends on Eigen and on the libigl 0.2.1 (http://igl.ethz.ch/projects/libigl/).

4) Matlab wrapper for the entire system

While the C++ code in 2) and 3) should be called directly from your C++ application to achieve good performances, we also provide two simple matlab wrappers that can be used to experiment with the method.

The two main methods are WA_forward.m and WA_inverse.m that solve the forward and inverse problem respectively. We also included a simple command-line test application WA_check.m and an interactive application that locally parametrize a mesh WA_starthere.m.

------ Dependencies

Matlab: Fast-marching toolbox (http://www.mathworks.com/matlabcentral/fileexchange/6110-toolbox-fast-marching), OpenMesh(http://www.openmesh.org).
C++: Eigen 3.2.0 (http://eigen.tuxfamily.org/), libigl(http://igl.ethz.ch/projects/libigl/)

All the dependencies (except OpenMesh) are included in the package for convenience.

Version 1.0 - August 14th 2013