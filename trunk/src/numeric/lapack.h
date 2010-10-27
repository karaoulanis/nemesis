/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]  *
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

#ifndef _LAPACK_H
#define _LAPACK_H

#include <acml.h>

#if defined(__GNUC__)
	#define dgesv	dgesv_
	#define dspsv	dspsv_
	#define dgbsv	dgbsv_
	#define dsyev	dsyev_
	#define dgetrf	dgetrf_
	#define dgetrs	dgetrs_
	#define dsptrf	dsptrf_
	#define dsptrs	dsptrs_
	#define dgbtrf	dgbtrf_
	#define dgbtrs	dgbtrs_
	#define dsygv	dsygv_
	#define dggev	dggev_
#endif

//#define __ACML__
#if defined(__ACML__)
	#define dgesv	DGESV
	#define dspsv	DSPSV
	#define dgbsv	DGBSV
	#define dsyev	DSYEV
	#define dgetrf	DGETRF
	#define dgetrs	DGETRS
	#define dsptrf	DSPTRF
	#define dsptrs	DSPTRS
	#define dgbtrf	DGBTRF
	#define dgbtrs	DGBTRS
	#define dsygv	DSYGV
	#define dggev	DGGEV
#endif

//extern "C" void dgesv(int *N, int *NRHS, double *A, int *LDA, 
//			      int *IPIV, double *B, int *LDB, int *INFO);
//extern "C" void dspsv(char* UPLO, int *N, int *NRHS, double *AP, 
//			      int *IPIV, double *B, int *LDB, int *INFO);
//extern "C" void dgbsv(int *N, int* KL, int* KU, int *NRHS, double *AB, int *LDAB, 
//			      int *IPIV, double *B, int *LDB, int *INFO);
//extern "C" void dsyev(char* JOBZ,char* UPLO,int* N,double* A,int* LDA,
//					  double* W,double* WORK,int* LWORK,int* INFO);
// 
// // Generic
//extern "C" void dgetrf(int* M,int* N,double* A,int* LDA,int* IPIV,int* INFO);
//extern "C" void dgetrs(char* TRANS,int* N,int* NRHS,double* A,int* LDA,
//					int* IPIV,double* B,int* LDB,int* INFO);
//
// 
////Symmetric
//extern "C" void dsptrf(char* UPLO,int* N,double* AP,int* IPIV,int* INFO);
//extern "C" void dsptrs(char* UPLO,int* N,int* NRHS,double* AP,
//						int* IPIV,double* B,int* LDB,int* INFO);
// 
//// Band
//extern "C" void dgbtrf(int* M,int* N,int* KL,int* KU,double* AB,
//					int* LDAB,int* IPIV,int* INFO);
//extern "C" void dgbtrs(char* TRANS,int* N,int* KL,int* KU,int* NRHS,double* AB,
//					int* LDAB,int* IPIV,double* B,int* LDB,int* INFO);
//
//// Generalized eigenvalue
//extern "C" void dsygv(int* ITYPE,char* JOBZ,char* UPLO,int* N,double* A,
//					  int* LDA,double* B,int* LDB,double* W,double* WORK,int* LWORK,int* INFO);
//extern "C" void dggev(char* JOBVL,char* JOBVR,int* N,double* A,
//					  int* LDA,double* B,int* LDB,
//					  double* ALPHAR,double* ALPHAI,double* BETA,
//					  double* VL,int* LDVL,double* VR,int* LDVR,double* WORK,int* LWORK,int* INFO);

#endif
