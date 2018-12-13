Fast 3x3 SVD
===========

This is an implementation of the method described in <a href="http://pages.cs.wisc.edu/~sifakis/papers/SVD_TR1690.pdf">"Computing the Singular Value Decomposition of 3x3 matrices with minimal branching and elementary floating point operations"</a>. I implemented this as part of <a href="http://wyegelwel.github.io/snow/">a group project</a> for a computer graphics course. 

Execution time per svd call on the CPU is about 2.0 microseconds. Tested on a AMD Phenom(tm) II X4 965 Processor. 

Execution time on the GPU is about 174 microseconds. Tested on a NVIDIA GeForce GTX 460 (profiled using nvvp).

Also included are routines for diagonalization / QR decomposition of 3x3 matrices, which may be useful in their own right. 


##Usage

Just include the header file and you are good to go! 

```C++

#include "svd3.h"
float a11, a12, a13, a21, a22, a23, a31, a32, a33;

a11= -0.558253; a12 = -0.0461681; a13 = -0.505735;
a21 = -0.411397; a22 = 0.0365854; a23 = 0.199707;
a31 = 0.285389; a32 =-0.313789; a33 = 0.200189;

float 	u11, u12, u13, 
		u21, u22, u23, 
		u31, u32, u33;

float 	s11, s12, s13, 
		s21, s22, s23, 
		s31, s32, s33;

float 	v11, v12, v13, 
		v21, v22, v23, 
		v31, v32, v33;

svd(a11, a12, a13, a21, a22, a23, a31, a32, a33,
    u11, u12, u13, u21, u22, u23, u31, u32, u33,
    s11, s12, s13, s21, s22, s23, s31, s32, s33,
    v11, v12, v13, v21, v22, v23, v31, v32, v33);

```

See the included Mathematica notebook for derivations of numerical shortcuts.

## License
MIT License, Eric V. Jang 2014