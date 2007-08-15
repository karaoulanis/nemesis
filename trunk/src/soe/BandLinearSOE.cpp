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

#include <BandLinearSOE.h>
#include <boost/numeric/ublas/io.hpp>


BandLinearSOE::BandLinearSOE()
	:SOE()
{
	myTag=TAG_SOE_LINEAR_BAND;
	lowerBandwidth=0;
	upperBandwidth=0;
	nRows=0;
}
BandLinearSOE::~BandLinearSOE()
{
}
int BandLinearSOE::
insertMatrixIntoA(const Matrix& Ke,const IDContainer& EFTable,double factor)
{
	isLUFactored=false;
	for(unsigned i=0;i<EFTable.size();i++)
		for(unsigned j=0;j<EFTable.size();j++)
		{
			if((EFTable[i])<0) continue;
			if((EFTable[j])<0) continue;
			A[EFTable[j]*nRows+(lowerBandwidth+upperBandwidth+EFTable[i]-EFTable[j])]+=factor*Ke(i,j);
		}
	return 0;
}
void BandLinearSOE::setTheSize()
{
	UndirectedGraph G;
	pA->getModel()->getUndirectedGraph(G);
	lowerBandwidth=bandwidth(G);
	upperBandwidth=lowerBandwidth;
		  
	// If the size has not changed do not resize arrays
	if(theSize==pA->getModel()->getnEquations()) return;
	else theSize=pA->getModel()->getnEquations();
	nRows=2*lowerBandwidth+upperBandwidth+1;
	A.resize(nRows*theSize);
	B.resize(theSize);
	X.resize(theSize);
	IPIV.resize(theSize);

	double d=(nRows+2)*theSize*sizeof(double);
	double i=theSize*sizeof(int);
	printf("soe: Allocated %.2fmb of memory for %ddofs.\n",(d+i)/(1024*1024),theSize);
}
void BandLinearSOE::print()
{
	for(int i=0;i<theSize;i++)
	{
		for(int j=0;j<theSize;j++)	
			cout<<A[j*nRows+(lowerBandwidth+upperBandwidth+i-j)]<<' ';
        cout<<" | "<<B[i]<<endl;
	}
}
int BandLinearSOE::getEigenSign()
{
	for(int i=0;i<theSize;i++) 
		if(fabs(A[lowerBandwidth+upperBandwidth+i*theSize])<1e-12)	return 0;
		else if(A[lowerBandwidth+upperBandwidth+i*theSize]<0)		return -1;
	return 1;
}
int BandLinearSOE::solve()
{
	// Input data
	int N = theSize;
	int KL=lowerBandwidth;
	int KU=upperBandwidth;
	int NRHS = 1;
	int INFO;
	int LDAB=2*KL+KU+1;
	int LDB=N;
	char c='N';
	X=B;

	// Compute the LU factorization of the band matrix A.
	if(!isLUFactored)
	{
		dgbtrf(&N,&N,&KL,&KU,&A[0],&LDAB,&IPIV[0],&INFO);
		if(INFO!=0)
			throw SException("[nemesis:%d] %s",1103,"SOE: lapack DGBTRF failed.\n");
		isLUFactored=true;
	}
	// Solve the system A*X = B, overwriting B with X.
	dgbtrs(&c,&N,&KL,&KU,&NRHS,&A[0],&LDAB,&IPIV[0],&X[0],&LDB,&INFO);
    return 0;
}
