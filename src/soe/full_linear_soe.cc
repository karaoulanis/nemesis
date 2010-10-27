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

#include "soe/FullLinearSOE.h"

FullLinearSOE::FullLinearSOE()
	:SOE()
{
	myTag=TAG_SOE_LINEAR_FULL;
}
FullLinearSOE::~FullLinearSOE()
{
}
int FullLinearSOE::
insertMatrixIntoA(const Matrix& Ke,const IDContainer& EFTable,double factor)
{
	isLUFactored=false;
	for(unsigned i=0;i<EFTable.size();i++)
		for(unsigned j=0;j<EFTable.size();j++)
		{
			if((EFTable[i])<0) continue;
			if((EFTable[j])<0) continue;
			A[EFTable[j]*theSize+EFTable[i]]+=factor*Ke(i,j);
		}
	return 0;
}
void FullLinearSOE::setTheSize()
{
	// If the size has not changed do not resize arrays
	if(theSize==pA->getModel()->getnEquations()) return;
	else theSize=pA->getModel()->getnEquations();
	A.resize(theSize*theSize);
	B.resize(theSize);
	X.resize(theSize);
	IPIV.resize(theSize);

	double d=((int)(theSize*theSize+2*theSize))*sizeof(double);
	double i=theSize*sizeof(int);
	printf("soe: Allocated %6.2fmb of memory for %ddofs.\n",(d+i)/(1024*1024),theSize);
}
void FullLinearSOE::print()
{
	for(int i=0;i<theSize;i++)
	{
		for(int j=0;j<theSize;j++)	cout<<A[j*theSize+i]<<' ';
        cout<<" | "<<B[i]<<endl;
	}
}
int FullLinearSOE::getEigenSign()
{
	for(int i=0;i<theSize;i++) 
		if(fabs(A[i*(theSize+1)])<1e-12)	return 0;
		else if(A[i*(theSize+1)]<0)			return -1;
	return 1;
}
int FullLinearSOE::solve()
{
	int N = theSize;
	int NRHS = 1;
	int LDA = N;
	int LDB = N;
	int INFO;
	char c='N';
//	this->print();
//	for(int i=0;i<theSize;i++)
//	{
//        cout<<B[i]<<endl;
//	}
//	cout<<endl;
	X=B;
	// Compute the LU factorization of the band matrix A.
	if(!isLUFactored)
	{
		dgetrf(&N,&N,&A[0],&LDA,&IPIV[0],&INFO);
		if(INFO!=0)	
			throw SException("[nemesis:%d] %s",1101,"SOE: lapack DGETRF failed.\n");
		isLUFactored=true;
	}
	// Solve the system A*X = B, overwriting B with X.
	dgetrs(&c,&N,&NRHS,&A[0],&LDA,&IPIV[0],&X[0],&LDB,&INFO,1);
	//for(int i=0;i<theSize;i++) cout<<B[i]<<endl; cout<<endl;
    return 0;
}
