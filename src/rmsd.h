/*
 *  *******************************************************************
 *
 *  rmsd.h 
 *  (c) 2005 Bosco K Ho
 * 
 *  Implementation of the Kabsch algorithm to find the RMSD, and 
 *  the least-squares rotation matrix for a superposition between 
 *  two sets of vectors.
 *
 *  This implementation is completely self-contained. No other dependencies.
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

#ifndef _rmsd_h_
#define _rmsd_h_

#ifdef __cplusplus
extern "C"
{
#endif

	/*
	 * calculate_rotation_rmsd()
	 *
	 *   given two lists of x,y,z coordinates, constructs
	 *    - mov_com: the centre of mass of the mov list
	 *    - mov_to_ref: vector between the com of mov and ref
	 *    - U: the rotation matrix for least-squares, usage of
	 *         of the matrix U[3][3] is
	 *           for (i=0; i<3; i++)
	 *           {
	 *             rotated_v[i] = 0.0;
	 *             for (j=0; j<3; j++)
	 *               rotated_v[i] += U[i][j] * v[j];
	 *           }
	 *    - rmsd: measures similarity between the vectors
	 */
	void calculate_rotation_rmsd(double ref_xlist[][3],
		double mov_xlist[][3],
		int n_list,
		double mov_com[3],
		double mov_to_ref[3],
		double U[3][3],
		double* rmsd);

	/*
	 * fast_rmsd()
	 *
	 * Fast calculation of rmsd w/o calculating a rotation matrix,
	 * adapted from the BTK by Chris Saunders 11/2002.
	 */
	void fast_rmsd(double ref_xlist[][3],
		double mov_xlist[][3],
		int n_list,
		double* rmsd);

#ifdef __cplusplus
}
#endif

#endif
