/*
 *  **************************************************************************
 *
 *  rmsd.c 
 *  (c) 2005 Bosco K Ho
 * 
 *  Implementation of the Kabsch algorithm to find the least-squares
 *  rotation matrix for a superposition between two sets of vectors.
 *  The jacobi transform of the diagonalization of a symmetric matrix
 *  is taken from Numerical Recipes. This piece of code is 
 *  self-contained and *does not require* another library.
 *
 *  **************************************************************************
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation; either version 2.1 of the License, or (at
 *  your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,  but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details. 
 *  
 *  You should have received a copy of the GNU Lesser General Public License 
 *  along with this program; if not, write to the Free Software Foundation, 
 *  Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 *  **************************************************************************
 * 
 */
 
#include <stdio.h>
#include <math.h>
#include "rmsd.h"



/* vector functions using c arrays */

void normalize(double a[3])
{
  double  b;

  b = sqrt((double)(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
  a[0] /= b;
  a[1] /= b;
  a[2] /= b;
}



double dot(double a[3], double b[3])
{
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}



static void cross(double a[3], double b[3], double c[3])
{
  a[0] = b[1]*c[2] - b[2]*c[1];
  a[1] = b[2]*c[0] - b[0]*c[2];
  a[2] = b[0]*c[1] - b[1]*c[0];
}



/*
 * setup_rotation() 
 *
 *      given two lists of x,y,z coordinates, constructs
 * the correlation R matrix and the E value needed to calculate the
 * least-squares rotation matrix.
 */
void setup_rotation(double ref_xlist[][3],
                    double mov_xlist[][3], 
                    int n_list,
                    double mov_com[3],
                    double mov_to_ref[3],
                    double R[3][3],
                    double* E0)
{
  int i, j, n;
  double ref_com[3];

  /* calculate the centre of mass */
  for (i=0; i<3; i++)
  { 
    mov_com[i] = 0.0;
    ref_com[i] = 0.0;
  }
  
  for (n=0; n<n_list; n++) 
    for (i=0; i<3; i++)
    { 
      mov_com[i] += mov_xlist[n][i];
      ref_com[i] += ref_xlist[n][i];
    }
    
  for (i=0; i<3; i++)
  {
    mov_com[i] /= n_list;
    ref_com[i] /= n_list;
    mov_to_ref[i] = ref_com[i] - mov_com[i];
  }

  /* shift mov_xlist and ref_xlist to centre of mass */
  for (n=0; n<n_list; n++) 
    for (i=0; i<3; i++)
    { 
      mov_xlist[n][i] -= mov_com[i];
      ref_xlist[n][i] -= ref_com[i];
    }

  /* initialize */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++) 
      R[i][j] = 0.0;
  *E0 = 0.0;

  for (n=0; n<n_list; n++) 
  {
    /* 
     * E0 = 1/2 * sum(over n): y(n)*y(n) + x(n)*x(n) 
     */
    for (i=0; i<3; i++)
      *E0 +=  mov_xlist[n][i] * mov_xlist[n][i]  
            + ref_xlist[n][i] * ref_xlist[n][i];
    
    /*
     * correlation matrix R:   
     *   R[i,j) = sum(over n): y(n,i) * x(n,j)  
     *   where x(n) and y(n) are two vector sets   
     */
    for (i=0; i<3; i++)
    {
      for (j=0; j<3; j++)
        R[i][j] += mov_xlist[n][i] * ref_xlist[n][j];
    }
  }
  *E0 *= 0.5;
  }



#define ROTATE(a,i,j,k,l) { g = a[i][j]; \
                            h = a[k][l]; \
                            a[i][j] = g-s*(h+g*tau); \
                            a[k][l] = h+s*(g-h*tau); }
/*   
 * jacobi3
 *
 *    computes eigenval and eigen_vec of a real 3x3
 * symmetric matrix. On output, elements of a that are above 
 * the diagonal are destroyed. d[1..3] returns the 
 * eigenval of a. v[1..3][1..3] is a matrix whose 
 * columns contain, on output, the normalized eigen_vec of
 * a. n_rot returns the number of Jacobi rotations that were required.
 */
int jacobi3(double a[3][3], double d[3], double v[3][3], int* n_rot)
{
  int count, k, i, j;
  double tresh, theta, tau, t, sum, s, h, g, c, b[3], z[3];

  /*Initialize v to the identity matrix.*/
  for (i=0; i<3; i++) 
  { 
    for (j=0; j<3; j++) 
      v[i][j] = 0.0;
    v[i][i] = 1.0;
  }

  /* Initialize b and d to the diagonal of a */
  for (i=0; i<3; i++) 
    b[i] = d[i] = a[i][i];

  /* z will accumulate terms */
  for (i=0; i<3; i++) 
    z[i] = 0.0; 
  
  *n_rot = 0;

  /* 50 tries */
  for (count=0; count<50; count++)     
  {

    /* sum off-diagonal elements */
    sum = 0.0;
    for (i=0; i<2; i++) 
    {
      for (j=i+1; j<3; j++)
         sum += fabs(a[i][j]);
    }

    /* if converged to machine underflow */
    if (sum == 0.0) 
      return(1);

    /* on 1st three sweeps... */
    if (count < 3) 
      tresh = sum * 0.2 / 9.0;    
    else       
      tresh = 0.0;      

    for (i=0; i<2; i++) 
    {
      for (j=i+1; j<3; j++) 
      {
        g = 100.0 * fabs(a[i][j]);

        /*  after four sweeps, skip the rotation if
         *   the off-diagonal element is small 
         */
        if ( count > 3  &&  fabs(d[i])+g == fabs(d[i])
              &&  fabs(d[j])+g == fabs(d[j]) ) 
        {
          a[i][j] = 0.0;
        } 
        else if (fabs(a[i][j]) > tresh) 
        {
          h = d[j] - d[i];
          
          if (fabs(h)+g == fabs(h))
          {
            t = a[i][j] / h;
          }
          else 
          {
            theta = 0.5 * h / (a[i][j]);
            t = 1.0 / ( fabs(theta) +
                        (double)sqrt(1.0 + theta*theta) );
            if (theta < 0.0) 
              t = -t;
          }
          
          c = 1.0 / (double) sqrt(1 + t*t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[i][j];

          z[i] -= h;
          z[j] += h;
          d[i] -= h;
          d[j] += h;

          a[i][j] = 0.0;

          for (k=0; k<=i-1; k++) 
            ROTATE(a, k, i, k, j)

          for (k=i+1; k<=j-1; k++) 
            ROTATE(a, i, k, k, j)

          for (k=j+1; k<3; k++) 
            ROTATE(a, i, k, j, k)

          for (k=0; k<3; k++) 
            ROTATE(v, k, i, k, j)

          ++(*n_rot);
        }
      }
    }

    for (i=0; i<3; i++) 
    {
      b[i] += z[i];
      d[i] = b[i];
      z[i] = 0.0;
    }
  }

  printf("Too many iterations in jacobi3\n");
  return (0);
}  



/* 
 * diagonalize_symmetric 
 *
 *    Diagonalize a 3x3 matrix & sort eigenval by size
 */
int diagonalize_symmetric(double matrix[3][3], 
                          double eigen_vec[3][3], 
                          double eigenval[3])
{
  int n_rot, i, j, k;
  double vec[3][3];
  double val; 
  
  if (!jacobi3(matrix, eigenval, vec, &n_rot)) 
  {
    printf("convergence failed\n");
    return (0);
  }

  /* sort solutions by eigenval */
  for (i=0; i<3; i++) 
  {
    k = i;
    val = eigenval[i];
    
    for (j=i+1; j<3; j++)
      if (eigenval[j] >= val)
      { 
        k = j;
        val = eigenval[k];
      }
       
    if (k != i) 
    {
      eigenval[k] = eigenval[i];
      eigenval[i] = val;
      for (j=0; j<3; j++) 
      {
        val = vec[j][i];
        vec[j][i] = vec[j][k];
        vec[j][k] = val;
      }
    }
  }

  /* transpose such that first index refers to solution index */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      eigen_vec[i][j] = vec[j][i];

  return (1);
}



/*
 * calculate_rotation_matrix() 
 *
 *   calculates the rotation matrix U and the
 * rmsd from the R matrix and E0:
 */
int calculate_rotation_matrix(double R[3][3],
                              double U[3][3], 
                              double E0,
                              double* residual)
{
  int i, j, k;
  double Rt[3][3], RtR[3][3];
  double left_eigenvec[3][3], right_eigenvec[3][3], eigenval[3];
  double v[3];
  double sigma;

  /* build Rt, transpose of R  */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      Rt[i][j] = R[j][i];

  /* make symmetric RtR = Rt X R */
  for (i=0; i<3; i++) 
    for (j=0; j<3; j++)
    {
      RtR[i][j] = 0.0;
      for (k = 0; k<3; k++)
        RtR[i][j] += Rt[k][i] * R[j][k];
    }

  if (!diagonalize_symmetric(RtR, right_eigenvec, eigenval))
    return(0);

  /* right_eigenvec's should be an orthogonal system but could be left
   * or right-handed. Let's force into right-handed system.
   */
  cross(&right_eigenvec[2][0], &right_eigenvec[0][0], &right_eigenvec[1][0]);

  /* From the Kabsch algorithm, the eigenvec's of RtR
   * are identical to the right_eigenvec's of R.
   * This means that left_eigenvec = R x right_eigenvec 
   */
  for (i=0; i<3; i++) 
    for (j=0; j<3; j++) 
      left_eigenvec[i][j] = dot(&right_eigenvec[i][0], &Rt[j][0]);

  for (i=0; i<3; i++) 
    normalize(&left_eigenvec[i][0]);

  /* 
   * Force left_eigenvec[2] to be orthogonal to the other vectors.
   * First check if the rotational matrices generated from the 
   * orthogonal eigenvectors are in a right-handed or left-handed
   * co-ordinate system - given by sigma. Sigma is needed to
   * resolve this ambiguity in calculating the RMSD.
   */
  cross(v, &left_eigenvec[0][0], &left_eigenvec[1][0]);
  if (dot(v, &left_eigenvec[2][0]) < 0.0)
    sigma = -1.0;
  else 
    sigma = 1.0;
  for (i=0; i<3; i++)
    left_eigenvec[2][i] = v[i]; 

  /* calc optimal rotation matrix U that minimises residual */
  for (i=0;i<3; i++)
    for (j=0; j<3; j++) 
    {
      U[i][j] = 0.0;
      for (k=0; k<3; k++)
        U[i][j] += left_eigenvec[k][i] * right_eigenvec[k][j];
    }
    
  *residual = E0 - (double) sqrt(fabs(eigenval[0])) 
                 - (double) sqrt(fabs(eigenval[1]))
                 - sigma * (double) sqrt(fabs(eigenval[2]));

  return (1);
}



void calculate_rotation_rmsd(double ref_xlist[][3],
                             double mov_xlist[][3], 
                             int n_list,
                             double mov_com[3],
                             double mov_to_ref[3],
                             double U[3][3],
                             double* rmsd)
{
  double Eo, residual;
  double R[3][3];
  
  setup_rotation(ref_xlist, mov_xlist, n_list, 
                 mov_com, mov_to_ref, R, &Eo);
  calculate_rotation_matrix(R, U, Eo, &residual);
  
  residual = fabs(residual); /* avoids the awkward case of -0.0 */
  *rmsd = sqrt( fabs((double) (residual)*2.0/((double)n_list)) ); 
}
 
 

/*
 * Fast calculation of rmsd w/o calculating a rotation matrix.
 *
 *   Chris Saunders 11/2002 - Fast rmsd calculation by the method of 
 * Kabsch 1978, where the required eigenvalues are found by an 
 * analytical, rather than iterative, method to save time. 
 * The cubic factorization used to accomplish this only produces 
 * stable eigenvalues for the transpose(R]*R matrix of a typical 
 * protein after the whole matrix has been normalized. Note that 
 * the normalization process used here is completely empirical 
 * and that, at the present time, there are **no checks** or 
 * warnings on the quality of the (potentially unstable) cubic 
 * factorization. 
 *
 */
#define PI 3.14159265358979323846
void fast_rmsd(double ref_xlist[][3],
               double mov_xlist[][3], 
               int n_list,
               double* rmsd)
{ 
  double R[3][3];
  double d0,d1,d2,e0,e1,f0;
  double omega;
  double mov_com[3];
  double mov_to_ref[3];

  /* cubic roots */
  double r1,r2,r3;
  double rlow;
  
  double v[3];
  double Eo, residual;
    
  setup_rotation(ref_xlist, mov_xlist, n_list, 
                 mov_com, mov_to_ref, R, &Eo);
  
  /* 
   * check if the determinant is greater than 0 to
   * see if R produces a right-handed or left-handed
   * co-ordinate system.
   */
  cross(v, &R[1][0], &R[2][0]);
  if (dot(&R[0][0], v) > 0.0)
    omega = 1.0;
  else
    omega = -1.0;

  /*
   * get elements we need from tran(R) x R 
   *  (funky matrix naming relic of first attempt using pivots)
   *          matrix = d0 e0 f0
   *                      d1 e1
   *                         d2
   * divide matrix by d0, so that cubic root algorithm can handle it 
   */
   
  d0 =  R[0][0]*R[0][0] + R[1][0]*R[1][0] + R[2][0]*R[2][0];

  d1 = (R[0][1]*R[0][1] + R[1][1]*R[1][1] + R[2][1]*R[2][1])/d0;
  d2 = (R[0][2]*R[0][2] + R[1][2]*R[1][2] + R[2][2]*R[2][2])/d0;

  e0 = (R[0][0]*R[0][1] + R[1][0]*R[1][1] + R[2][0]*R[2][1])/d0;
  e1 = (R[0][1]*R[0][2] + R[1][1]*R[1][2] + R[2][1]*R[2][2])/d0;

  f0 = (R[0][0]*R[0][2] + R[1][0]*R[1][2] + R[2][0]*R[2][2])/d0;

  /* cubic roots */
  {
    double B, C, D, q, q3, r, theta;
    /*
     * solving for eigenvalues as det(A-I*lambda) = 0
     * yeilds the values below corresponding to:
     * lambda**3 + B*lambda**2 + C*lambda + D = 0
     *   (given that d0=1.)
     */
    B = -1.0 - d1 - d2;
    C = d1 + d2 + d1*d2 - e0*e0 - f0*f0 - e1*e1;
    D = e0*e0*d2 + e1*e1 + f0*f0*d1 - d1*d2 - 2*e0*f0*e1;

    /* cubic root method of Viete with all safety belts off */
    q = (B*B - 3.0*C) / 9.0;
    q3 = q*q*q;
    r = (2.0*B*B*B - 9.0*B*C + 27.0*D) / 54.0;
    theta = acos(r/sqrt(q3));
    r1 = r2 = r3 = -2.0*sqrt(q);
    r1 *= cos(theta/3.0);
    r2 *= cos((theta + 2.0*PI) / 3.0);
    r3 *= cos((theta - 2.0*PI) / 3.0);
    r1 -= B / 3.0;
    r2 -= B / 3.0;
    r3 -= B / 3.0;
  }

  /* undo the d0 norm to get eigenvalues */
  r1 = r1*d0;
  r2 = r2*d0;
  r3 = r3*d0;

  /* set rlow to lowest eigenval; set other two to r1,r2 */
  if (r3<r1 && r3<r2)
  {
    rlow = r3;
  }
  else if (r2<r1 && r2<r3)
  {
    rlow = r2; 
    r2 = r3;
  } 
  else 
  { 
    rlow = r1; 
    r1 = r3;
  }

  residual = Eo - sqrt(r1) - sqrt(r2) - omega*sqrt(rlow);
  *rmsd = sqrt( (double) residual*2.0 / ((double) n_list) ); 
}



