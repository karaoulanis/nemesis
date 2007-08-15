/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 3, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License         *
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//*****************************************************************************

#ifndef _LU_H
#define _LU_H

namespace LU
{
	///@todo check and comment the algorithm
	inline int decomposition(double* data,int n,double* vv,int* index,double& d) 
	{
		int i,imax,j,k;
		double big,sum,temp;
		d=1;
		//search for the largest element in each row; save the scaling in the 
		//temporary array vv and return -1 if the matrix is singular
		for(i=0;i<n;i++) 
		{
			big=0.;
			for(j=0;j<n;j++) if((temp=fabs(data[i*n+j]))>big) big=temp;
			if(big==0.) return -1;
			vv[i]=big;
		}
		/* the main loop for the Crout's algorithm */
		for(j=0;j<n;j++) 
		{
			/* this is the part a) of the algorithm except for i==j */
			for(i=0;i<j;i++) 
			{
				sum=data[i*n+j];
				for(k=0;k<i;k++) sum-=data[i*n+k]*data[k*n+j];
				data[i*n+j]=sum;
			}
			/* initialize for the search for the largest pivot element */
			big=0.;
			imax=j;
			/* this is the part a) for i==j and part b) for i>j + pivot search */
			for(i=j;i<n;i++) 
			{
				sum=data[i*n+j];
				for(k=0;k<j;k++) sum-=data[i*n+k]*data[k*n+j];
				data[i*n+j]=sum;
				/* is the figure of merit for the pivot better than the best so far? */
				if((temp=vv[i]*fabs(sum))>=big) 
				{
					big=temp;
					imax=i;
				}
			}
			/* interchange n, if needed, change parity and the scale factor */
			if(imax!=j) 
			{
				for(k=0;k<n;k++) 
				{
					temp=data[imax*n+k];
					data[imax*n+k]=data[j*n+k];
					data[j*n+k]=temp;
				}
				d*=-1.0;
				vv[imax]=vv[j];
			}
			/* store the index */
			index[j]=imax;
			/* if the pivot element is zero, the matrix is singular but for some 
			applications a tiny number is desirable instead */
			if(data[j*n+j]==0.) data[j*n+j]=num::eps;
			/* finally, divide by the pivot element */
			if(j<n-1) 
			{
				temp=1./data[j*n+j];
				for(i=j+1;i<n;i++) data[i*n+j]*=temp;
			}
		}
		return 0;
	}
	inline void backsubstitution(double* data,int n,int* index,double* b) 
	{
		int i,j,ip,ii=-1;
		double sum;
		// First step of backsubstitution; the only wrinkle is to unscramble 
		// the permutation order. Note: the algorithm is optimized for a 
		// possibility of large amount of zeroes in b
		for(i=0;i<n;i++) 
		{
			ip=index[i];
			sum=b[ip];
			b[ip]=b[i];
			if(ii>=0) for(j=ii;j<i;j++) sum-=data[i*n+j]*b[j];
			else if(sum) ii=i;	// a nonzero element encounted
			b[i]=sum;
		}
		// the second step
		for(i=n-1;i>=0;i--) 
		{
			sum=b[i];
			for(j=i+1;j<n;j++) sum-=data[i*n+j]*b[j];
			b[i]=sum/data[i*n+i];
		}
	}
	inline void inverse(double* data,int n,int* index,double* inv,double *col)
	{
		for(int j=0;j<n;j++) 
		{
			for(int i=0;i<n;i++) col[i]=0.;
			col[j]=1.;
			backsubstitution(data,n,index,col);
			for(int i=0;i<n;i++) inv[i*n+j]=col[i];
		}
	}
	inline double determinant(double* data,int n,double& d) 
	{
		double det=d;
		for(int j=0;j<n;j++) det*=data[j*n+j];
		return det;
	}
}
#endif
