#include <iostream>
#include <glm/geometric.hpp>
#include <glm/gtc/epsilon.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/fast_square_root.hpp>
#include <glm/gtx/fast_trigonometry.hpp>
#include <glm/gtc/matrix_access.hpp>
#include "math.h"

#define _sqrtHalf 0.70710678
#define _gamma 5.828427124 // FOUR_GAMMA_SQUARED = sqrt(8)+3;
#define _cstar 0.923879532 // cos(pi/8)
#define _sstar 0.3826834323 // sin(p/8)

#define _EPSILON_ 1e-6
#define EPSILON _EPSILON_

#define EQ(a, b) (fabs((a) - (b)) < _EPSILON_)
#define NEQ(a, b) (fabs((a) - (b)) > _EPSILON_)

void printMat3(glm::mat3 mat) {
    // prints by rows
    for (int j=0; j<3; ++j) // g3d stores column-major
    {
        for (int i=0; i<3; ++i)
        {
            printf("%f   ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printQuat(glm::quat q)
{
    std::cout << q.w << "  " << q.x << "  " << q.y  << "  " << q.z <<std::endl;
    //std::cout << q[3] << "  " << q[0] << "  " << q[1]  << "  " << q[2] <<std::endl;
}

inline float accurateRSQRT(float x)
{
    // used in step 3
    // cool result from Lomont 2003
    float y = glm::fastSqrt(x);
    return y * (3-x*y*y)/2;
}

inline float accurateSQRT(float x) { return x * accurateRSQRT(x); } // used in step 3

void condSwap(bool c, float &X, float &Y)
{
    // used in step 2
    float Z = X;
    X = c ? Y : X;
    Y = c ? Z : Y;
}

// swapping functions for entire rows
void condSwap(bool c, glm::vec3 &X, glm::vec3 &Y)
{
    // used in step 2
    glm::vec3 Z = X;
    X = c ? Y : X;
    Y = c ? Z : Y;
}

void condNegSwap(bool c, glm::vec3 &X, glm::vec3 &Y)
{
    // used in step 2 and 3
    glm::vec3 Z = -X;
    X = c ? Y : X;
    Y = c ? Z : Y;
}

void condNegSwap(bool c, int c1, int c2, glm::quat &qV)
{
    // condNegSwap can be modified to operate well
    // on quaternion representation of V
    glm::quat qR;
    // qR = (1, 0, 0, c) for 1,2
    // what are the other quaternions corresponding to the other rot matrices?
    qV = qV * qR;
}

//void approximateGivensAngles(float a11, float a12, float a22, float &c, float &s)
//{
//    // this function is unused
//    // problem, second round the angles are flipped?
//    /*
//     * approximates givens angles for 2x2 submatrix. used in step 1
//     *
//     * A = (a11 a12
//     *      a12 a22)
//     *
//     * where 1,2 correspond to p,q values for givens rotation of a 3x3 matrix
//     */
//    float A2 = a12 * a12;
//    float B2 = (a11-a22)*(a11-a22);
//    bool b = A2 < B2;
//    //float w = rsqrt(A2+B2);
//    float w = glm::fastInverseSqrt(A2+B2);
//    s = b ? w*a12 : _sqrtHalf;
//    c = b ? w*(a11-a22) : _sqrtHalf;
//    // we can represent this rotation with quaternion
//    // (cos(theta/2),0,0,sin(theta/2)) for p,q = 1,2
//    // up to the user to build the appropriate quaternion or rotation matrix
//}

void approximateGivensQuaternion(float a11, float a12, float a22, float &ch, float &sh)
{
    /*
     * Given givens angle computed by approximateGivensAngles,
     * compute the corresponding rotation quaternion.
     */
    // used during second half of each jacobi iteration step?
    ch = 2*(a11-a22);
    sh = a12;
    bool b = _gamma*sh*sh < ch*ch;
    float w = glm::fastInverseSqrt(ch*ch+sh*sh);
    ch=b?w*ch:_cstar;
    sh=b?w*sh:_sstar;
    //qV=glm::quat(ch,0,0,sh)
    //up to the user to build the appropriate quaternion or unscaled rotation matrix
}

void jacobiConjugation(int p, int q, glm::mat3 &S, glm::quat &qV)
{
    //std::cout << "matrix before conjugation" << std::endl;
    //printMat3(S);
    // eliminate off-diagonal entries Spq, Sqp
    float ch,sh;
    approximateGivensQuaternion(S[0][0],S[1][0],S[1][1],ch,sh);

    // build rotation matrix Q
    glm::mat3 Q;
    float scale = ch*ch+sh*sh;
    float a = (ch*ch-sh*sh)/scale;
    float b = (2*sh*ch)/scale;
    
//    printf("a=%f, b=%f \n", a, b );

    Q[0][0] = a;    Q[1][0] = -b;
    Q[0][1] = b;    Q[1][1] = a;

    // perform the conjugation to annihilate S = Q' S Q
    S =glm::transpose(Q) * S * Q;


    // std::cout << "result of conjugation" << std::endl;
    // printMat3(S);
    // update cumulative rotation quaternion qV
    // tmp1, tmp2, tmp3
    glm::vec3 tmp(qV.x,qV.y,qV.z);
    tmp *= sh;
    sh *= qV.w;
    // printf("sh=%f\n",sh);
    // original
    qV *= ch;

    // this arranges such that for (p,q) = ((0,1),(1,2),(0,2)),
    // n = (0,1,2)
    int n = 2*q-p-2; // worst hack in the world

    int x = n;
    int y = (n+1)%3;
    int z = (n+2)%3;

    qV[z] += sh;
    qV.w -= tmp[z];
    qV[x] += tmp[y];
    qV[y] -= tmp[x];

    // re-arrange matrix for next iteration
    S = glm::mat3(S[1][1], S[1][2], S[1][0],
                  S[1][2], S[2][2], S[0][2],
                  S[1][0], S[0][2], S[0][0]);
}

void sortSingularValues(glm::mat3 &B, glm::mat3 &V)
{
    // used in step 2
    glm::vec3 b1 = glm::column(B,0); glm::vec3 v1 = glm::column(V,0);
    glm::vec3 b2 = glm::column(B,1); glm::vec3 v2 = glm::column(V,1);
    glm::vec3 b3 = glm::column(B,2); glm::vec3 v3 = glm::column(V,2);
    float rho1 = glm::length2(b1);
    float rho2 = glm::length2(b2);
    float rho3 = glm::length2(b3);
    bool c;

    c = rho1 < rho2;
    condNegSwap(c,b1,b2); condNegSwap(c,v1,v2);
    condSwap(c,rho1,rho2);

    c = rho1 < rho3;
    condNegSwap(c,b1,b3); condNegSwap(c,v1,v3);
    condSwap(c,rho1,rho3);

    c = rho2 < rho3;
    condNegSwap(c,b2,b3); condNegSwap(c,v2,v3);

    // re-build B,V
    B = glm::mat3(b1,b2,b3);
    V = glm::mat3(v1,v2,v3);
}


void QRGivensQuaternion(float a1, float a2, float &ch, float &sh)
{
    // a1 = pivot point on diagonal
    // a2 = lower triangular entry we want to annihilate
    float epsilon = EPSILON;

    float rho = sqrt(a1*a1 + a2*a2);

    sh = rho > epsilon ? a2 : 0;
    ch = fabs(a1) + fmax(rho,epsilon);
    bool b = a1 < 0;
    condSwap(b,sh,ch);
    float w = glm::inversesqrt(ch*ch+sh*sh);
    //float w = glm::fastInverseSqrt(ch*ch+sh*sh);
    ch *= w;
    sh *= w;
}

void QRDecomposition(glm::mat3 B, glm::mat3 &Q, glm::mat3 &R)
{
    // QR decomposition of 3x3 matrices using Givens rotations to
    // eliminate elements B21, B31, B32

    glm::quat qQ; // cumulative rotation
    glm::quat qU; // each Givens rotation in quaternion form
    glm::mat3 U;
    float ch, sh;
    QRGivensQuaternion(B[0][0],B[0][1],ch,sh);
    qU = glm::quat(ch,0,0,sh);
    U = glm::toMat3(qU);
    B = glm::transpose(U) * B;

    // update cumulative rotation
    qQ *= qU;

    // second givens rotation
    QRGivensQuaternion(B[0][0],B[0][2],ch,sh);
    qU = glm::quat(ch,0,-sh,0);
    U = glm::toMat3(qU);
    B = glm::transpose(U) * B;
    qQ *= qU;

    // third Givens rotation
    QRGivensQuaternion(B[1][1],B[1][2],ch,sh);
    qU = glm::quat(ch,sh,0,0);
    U = glm::toMat3(qU);
    B = glm::transpose(U) * B;
    qQ *= qU;
    // B has been transformed into R
    R = B;
    // qQ now contains final rotation for Q
    Q = glm::toMat3(qQ);
}

void jacobiEigenanalysis(glm::mat3 &S, glm::quat &qV)
{
    // wrapper function for the first step
    // solve symmetric eigenproblem using jacobi iteration
    // given a symmetric matrix S, diagonalize it
    // also returns the cumulative rotation as a quaternion
    qV = glm::quat(1,0,0,0);
    float scale;
    int p; int q;
    for(int sweep=0;sweep<4;sweep++)
    {
        // we wish to eliminate the maximum off-diagonal element
        // on every iteration, but cycling over all 3 possible rotations
        // in fixed order (p,q) = (1,2) , (2,3), (1,3) still has
        //  asymptotic convergence
        jacobiConjugation(0,1,S,qV);
        jacobiConjugation(1,2,S,qV);
        jacobiConjugation(0,2,S,qV);
    }
}


int main()
{

    // run a simpe test
    glm::mat3 A;
    // GLM stores matrices in column-major order so this initialization is the transpose of what we want
    A = glm::mat3(   -0.558253,    -0.0461681,    -0.505735,
                   -0.411397  ,   0.0365854  ,   0.199707,
                    0.285389  ,   -0.313789  ,   0.200189);

    A = glm::transpose(A);

    std::cout << "original matrix" << std::endl;
    printMat3(A);
/// 2. Symmetric Eigenanlysis
    // normal equations matrix
    glm::mat3 S = glm::transpose(A) * A;
    // std::cout << "normal equations matrix" << std::endl;
    // printMat3(S);
    glm::quat qV;
    jacobiEigenanalysis(S,qV);
    std::cout << "final S (diagonalized) " << std::endl;
    printMat3(S);
    std::cout << "cumulative rotation quaternion" << std::endl;
    printQuat(qV);


    glm::mat3 V = glm::toMat3(qV); // normalize qV, convert it to matrix
    //std::cout << "final rot matrix V" << std::endl;
    //printMat3(V);
    glm::mat3 B = A * V; // right-multiply A with V => left multiply for column major
    std::cout << "B=AV" << std::endl;
    //printMat3(B);
///// 3. Sorting the singular values (find V)
    sortSingularValues(B,V);
    //std::cout << "sorted B=AV" << std::endl;
    //printMat3(B);
    std::cout << "sorted V" << std::endl;
    printMat3(V);

//    // if columns 2-3 swapped, also update quaternion representation (i dont think we need to though since we're not tracking quats)
///// 4. QR factorization using Givens rotations (find U,S from B=AV)
    glm::mat3 U;
    glm::mat3 Sigma;
    QRDecomposition(B,U,Sigma);
    std::cout << "U" << std::endl;
    printMat3(U);
    std::cout << "Sigma" << std::endl;
    printMat3(Sigma);
    glm::mat3 USV = U * Sigma * glm::transpose(V);
    std::cout << "product USV'"<< std::endl;
    printMat3(USV);

    return 0;
}

