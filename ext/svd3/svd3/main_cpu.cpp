#include <stdio.h>
#include "svd3.h"
#include <time.h>

// some printing utilities
inline void printMat3(float a11, float a12, float a13,
               float a21, float a22, float a23,
               float a31, float a32, float a33)
{
    printf("%f %f %f \n", a11, a12, a13);
    printf("%f %f %f \n", a21, a22, a23);
    printf("%f %f %f \n", a31, a32, a33);
}

inline void printQuat(float * q)
{
    // print w,x,y,z
    printf("%f %f %f %f\n",q[3],q[0],q[1],q[2]);
}


int main(void)
{
    float a11, a12, a13, a21, a22, a23, a31, a32, a33;

    a11= -0.558253; a12 = -0.0461681; a13 = -0.505735;
    a21 = -0.411397; a22 = 0.0365854; a23 = 0.199707;
    a31 = 0.285389; a32 =-0.313789; a33 = 0.200189;

    // printf("Original Matrix:\n");
    // printMat3(a11, a12, a13, a21, a22, a23, a31, a32, a33);

    float u11, u12, u13, u21, u22, u23, u31, u32, u33;
    float s11, s12, s13, s21, s22, s23, s31, s32, s33;
    float v11, v12, v13, v21, v22, v23, v31, v32, v33;

    clock_t start, end;
    start = clock();
    for (int i=0; i<1e6; i++)
    {
        svd(a11, a12, a13, a21, a22, a23, a31, a32, a33,
        u11, u12, u13, u21, u22, u23, u31, u32, u33,
        s11, s12, s13, s21, s22, s23, s31, s32, s33,
        v11, v12, v13, v21, v22, v23, v31, v32, v33);
    }    
    end = clock();

    printf("Average SVD takes %f microseconds \n ", 1e6*(double(end - start) / 1e6 / CLOCKS_PER_SEC ) );
    
    printf("U:\n");
    printMat3(u11, u12, u13, u21, u22, u23, u31, u32, u33);
    printf("S:\n");
    printMat3(s11, s12, s13, s21, s22, s23, s31, s32, s33);
    printf("V:\n");
    printMat3(v11, v12, v13, v21, v22, v23, v31, v32, v33);

    float t11, t12, t13, t21, t22, t23, t31, t32, t33;
    multAB(u11, u12, u13, u21, u22, u23, u31, u32, u33,
           s11, s12, s13, s21, s22, s23, s31, s32, s33,
           t11, t12, t13, t21, t22, t23, t31, t32, t33);

    float m11, m12, m13, m21, m22, m23, m31, m32, m33;
    multAB(t11, t12, t13, t21, t22, t23, t31, t32, t33,
           v11, v21, v31, v12, v22, v32, v13, v23, v33,
           m11, m12, m13, m21, m22, m23, m31, m32, m33);

    printf("USV* : \n");
    printMat3(m11, m12, m13, m21, m22, m23, m31, m32, m33);

    return 0;

}

