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
//*****************************************************************************

#include <BFGS.h>

BFGS::BFGS(int m_)
:m(m_),isLineSearchActive(false)
{
	myTag=TAG_NONE;
	s.resize(m);
	y.resize(m);
}
BFGS::BFGS(int m_,double etaMin_,double etaMax_,double rTol_,int maxIter_)
:m(m_),etaMin(etaMin_),etaMax(etaMax_),rTol(rTol_),maxIter(maxIter_),isLineSearchActive(true)
{
	myTag=TAG_NONE;
	s.resize(m);
	y.resize(m);
}
BFGS::~BFGS()
{
}
int BFGS::solveStep(int n)
{
	int size=pA->getModel()->getnEquations();
	Vector resOld(size);
	Vector resNew(size);
	Vector q(size);
	Vector du(size);
	Vector a(m);
	Vector r(m);
	for(int j=0;j<m;j++)
	{
		s[j].resize(size,0.);
		y[j].resize(size,0.);
	}

	// Predictor phase
	pA->getControl()->formTangent();
	pA->getControl()->predict();
	pA->getConvergenceNorm()->newStep();
	du=pA->getSOE()->getX();
	resOld=pA->getSOE()->getB();
	resOld*=-1.;
	pA->getControl()->formResidual(pA->getControl()->getLambda());
	
	// Corrector phase
	int k=0;
	int check;
	while((check=pA->getConvergenceNorm()->update())>0)
	{
		resNew=pA->getSOE()->getB();
		resNew*=-1.;
		y[k]=resNew-resOld;
		s[k]=du;

		q=-resNew;
		//for(int i=k-1;i>=k-m;i--)				// for l-bfgs
		for(int i=k-1;i>=0;i--)
		{
			//cout<<"A "<<i<<" "<<k-i-1<<endl;	// for l-bfgs
			//if(i<=0) continue;				// for l-bfgs
			r[i+1]=1/(y[i+1]*s[i+1]);
			a[i]=r[i+1]*s[i+1]*q;
			q-=a[i]*y[i+1];
		}
		pA->getSOE()->setB(q);
		pA->getSOE()->solve();
		du=pA->getSOE()->getX();
		//for(int i=k-m;i<=k-1;i++)				// for l-bfgs
		for(int i=1;i<=k;i++)
		{
			//cout<<"B "<<i<<" "<<k-i-1<<endl;	// for l-bfgs
			//if(i<=0) continue;				// for l-bfgs
			double b=r[i]*y[i]*du;
			du+=s[i]*(a[i-1]-b);
		}

		pA->getSOE()->setX(du);
		pA->getControl()->correct();
		pA->getControl()->formResidual(pA->getControl()->getLambda());
		
		double s0=-s[k]*resNew;
		double s1=-s[k]*(pA->getSOE()->getB());
		if(isLineSearchActive) this->lineSearch(s0,s1,du);
		
		resOld=resNew;
		++k;
 	}
	return check;      
}
void BFGS::lineSearch(double s0,double sj,const Vector& du)
{
	Vector dx=du;
	int k=0;
	double eta=1.0;
	double etaP=1.0;
	
	double r0;
	if(!num::tiny(s0))	r0=fabs(sj/s0);
	else				r0=0.;
	double rj=r0;
	while(k<maxIter && rj>rTol)
	{
		eta=etaP*s0/(s0-sj);
		if(eta<etaMin) eta=etaMin;			
		if(eta>etaMax) eta=etaMax;
		if(rj >r0    ) eta=1.0;
		dx*=eta;
		pA->getSOE()->setX(dx);
		pA->getControl()->correct();
		pA->getControl()->formResidual(pA->getControl()->getLambda());
		sj=du*(pA->getSOE()->getB());
		rj=fabs(sj/s0);
		etaP=eta;
		++k;
	}
	dx=eta*du;
	pA->getSOE()->setX(dx);
 }
