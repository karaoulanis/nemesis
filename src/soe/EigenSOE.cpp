/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 2, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License along   *
*   with this program; if not, write to the Free Software Foundation, Inc.,   *
*   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.               *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//            C.G. Panagiotopoulos (pchr@civil.auth.gr)
//*****************************************************************************

#include <EigenSOE.h>

EigenSOE::EigenSOE()
	:SOE()
{
	myTag=TAG_NONE;
}
EigenSOE::~EigenSOE()
{
}
int EigenSOE::
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
int EigenSOE::
insertMatrixIntoM(const Matrix& Ke,const IDContainer& EFTable,double factor)
{
	isLUFactored=false;
	for(unsigned i=0;i<EFTable.size();i++)
		for(unsigned j=0;j<EFTable.size();j++)
		{
			if((EFTable[i])<0) continue;
			if((EFTable[j])<0) continue;
			M[EFTable[j]*theSize+EFTable[i]]+=factor*Ke(i,j);
		}
	return 0;
}
void EigenSOE::setTheSize()
{
	// If the size has not changed do not resize arrays
	if(theSize==pA->getModel()->getnEquations()) return;
	else theSize=pA->getModel()->getnEquations();
	A.resize(theSize*theSize); 
	M.resize(theSize*theSize);
	X.resize(theSize);
	ALPHAR.resize(theSize); 
	ALPHAI.resize(theSize); 
	BETA.resize(theSize); 
	VL.resize(theSize*theSize);
	VR.resize(theSize*theSize); 
	WORK.resize(8*theSize);

	double d=(4*theSize*theSize+9*theSize)*sizeof(double);
	double i=theSize*sizeof(int);
	printf("soe: Allocated %6.2fmb of memory.\n",(d+i)/(1024*1024));
}
void EigenSOE::print()
{
	cout<<"A"<<endl;
	for(int i=0;i<theSize;i++)
	{
		for(int j=0;j<theSize;j++)	cout<<A[j*theSize+i]<<' ';
		cout<<endl;
	}
	cout<<"M"<<endl;
	for(int i=0;i<theSize;i++)
	{
		for(int j=0;j<theSize;j++)	cout<<M[j*theSize+i]<<' ';
		cout<<endl;
	}
}
int EigenSOE::getEigenSign()
{
	return 1;
}
void EigenSOE::zeroM()
{
	for(unsigned i=0;i<M.size();i++) M[i]=0;
}
int EigenSOE::solve()
{
	char JOBVL='N';
	char JOBVR='V';
	int N = theSize;
	int LDA = N;
	int LDB = N;
	int LDVL = N;
	int LDVR = N;
	int LWORK = 8*N;
	int INFO;
	dggev(&JOBVL,&JOBVR,&N,&A[0],&LDA,&M[0],&LDB,&ALPHAR[0],&ALPHAI[0],&BETA[0],
					  &VL[0],&LDVL,&VR[0],&LDVR,&WORK[0],&LWORK,&INFO);
	if(INFO!=0)	
		throw SException("[nemesis:%d] %s",1110,"SOE: lapack DSYGV failed.");
	for(int i=0;i<theSize;i++) X[i]=ALPHAR[i]/BETA[i];
    return 0;
}
