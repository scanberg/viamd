#include <cuda.h>
#include "svd3_cuda.h"
#include <stdio.h>

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
  # error printf is only supported on devices of compute capability 2.0 and higher, please compile with -arch=sm_20 or higher
#endif

extern "C" void runCudaPart();

__host__ __device__ __forceinline__
void printMat3(float a11, float a12, float a13,
               float a21, float a22, float a23,
               float a31, float a32, float a33)
{
    printf("%f %f %f \n", a11, a12, a13);
    printf("%f %f %f \n", a21, a22, a23);
    printf("%f %f %f \n", a31, a32, a33);
}

__global__ void svd3_test()
{
    int tid = blockIdx.x;

    float a11, a12, a13, a21, a22, a23, a31, a32, a33;

    a11= -0.558253; a12 = -0.0461681; a13 = -0.505735;
    a21 = -0.411397; a22 = 0.0365854; a23 = 0.199707;
    a31 = 0.285389; a32 =-0.313789; a33 = 0.200189;

    float u11, u12, u13, u21, u22, u23, u31, u32, u33;
    float s11, s12, s13, s21, s22, s23, s31, s32, s33;
    float v11, v12, v13, v21, v22, v23, v31, v32, v33;

    svd(a11, a12, a13, a21, a22, a23, a31, a32, a33,
    u11, u12, u13, u21, u22, u23, u31, u32, u33,
    s11, s12, s13, s21, s22, s23, s31, s32, s33,
    v11, v12, v13, v21, v22, v23, v31, v32, v33);
}

void runCudaPart()
{
    // all your cuda code here
    svd3_test<<<1,1>>>(); // 5 blocks, 1 GPU thread each
    cudaDeviceSynchronize();
}
