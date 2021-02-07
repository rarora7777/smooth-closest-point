Code to compute offset-surfaces around a manifold triangle mesh, tetrahedral mesh the space between the inner and outer offset surfaces (the _shell_), and embed the tetrahedral mesh in a higher-dimensional space. This higher-dimensional embedding can then be used for smooth closest-point queries for object within the shell.

Based on [Panozzo et al. 2013], see below.

# Build Instructions

Dependencies: [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/download/) (7.1 or higher), [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page), [libigl](https://github.com/libigl/libigl/). [TetWild](https://github.com/Yixin-Hu/TetWild) is also recommended, but is not required for building the application.

Please modify the paths below based on your machine. `LIBIGL_DIR` should point to the libigl root folder and `Eigen3_DIR` should point to the folder containing `Eigen3Config.cmake`.

## Windows

```
mkdir build
cd build
cmake .. -DLIBIGL_DIR="D:/dev/libigl" -DEigen3_DIR="D:/dev/eigen3/build"
cmake --build . --config Release
```

## Linux/MacOS

```
mkdir build
cd build
cmake .. -DLIBIGL_DIR="/usr/local/libigl" -DEigen3_DIR="/usr/local/eigen3/build" -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

# Usage


## Pre-processing
For pre-processing (compute and embed shell tet mesh), the main point of entry is `computeAndSaveOffsetSurfaces` followed by `computeAndSaveEmbedding`.

- Download a manifold triangle mesh in Wavefront OBJ format to `./data/`. Let's call it `bunny.obj`.
- In MATLAB, run `computeAndSaveEmbedding(bunny)`.
- [Optional, but recommended] Use TetWild to improve mesh quality of the offset surfaces `./data/bunny_out_cgal.obj` and `./data/bunny_in_cgal.obj`. For objects that are thin everywhere, the inner offset surface may be absent. Please ensure that the TetWild-processed surfaces are located at `./data/bunny_out.obj` and `./data/bunny_in.obj` (if inner surface exists).
- In MATLAB, run `computeAndSaveEmbedding(bunny)`. This may take a while. The triangle and tetrahedral meshes are created in a simple text format at `./data/bunny_tri.txt` and `./data/bunny_tet.txt`. Use `readMesh('./data/bunny_tri.txt', 3, 8)` and `readMesh('./data/bunny_tet.txt', 4, 8)` to read and view these meshes, if required.


## Runtime

Call functions in `Phong.dll` (`Phong.so`/`Phong.dylib` on Linux/MacOS). Typical usage:

- Create a `Phong` object using `createPhongObject()`.
- For every point that needs to be projected, call `project()` or `projectBruteForce()`.
- Delete the `Phong` object using `deletePhongObject()`.

Please see `src/phong/Phong.h` for details. The exported functions are at the bottom of the file in an `extern "C"` block.


# License

In accordance with the license associated with Panozzo et al.'s original repo, this code is free to use for non-commercial applications, including my changes and additions.


# Stuff from Panozzo et al.'s original repo
------ Introduction

Reference implementation for the paper Weighted Averages on Surfaces
Daniele Panozzo, Ilya Baran, Olga Diamanti, Olga Sorkine
SIGGRAPH 2013

This code is free to use for non-commercial applications, please cite our paper if you use it in your project.

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
