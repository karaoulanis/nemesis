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

#include "material/drucker_prager_new3.h"
#include "material/drucker_prager_ys.h"
#include "material/tension_cutoff_ys.h"

Matrix DruckerPragerNew3::C(6,6,0.);
Matrix DruckerPragerNew3::C3(3,3,0.);

DruckerPragerNew3::DruckerPragerNew3()
{
}
DruckerPragerNew3::DruckerPragerNew3(int ID,int elasticID,double c,double phi,double psi,double Kc,double Kphi,double T)
:MultiaxialMaterial(ID,0.,0.)
{
	// Get the elastic part
	Material* p=pD->get<Material>(pD->getMaterials(),elasticID);
	if(p->getTag()!=TAG_MATERIAL_MULTIAXIAL_ELASTIC)
		throw SException("[nemesis:%d] %s",9999,"Multiaxial elastic material expected.");
	myElastic=static_cast<MultiaxialMaterial*>(p)->getClone();
	MatParams[30]=myElastic->getParam(30);
	MatParams[31]=myElastic->getParam(31);

	// Material properties
	MatParams[0]=c;
	MatParams[1]=phi;
	MatParams[2]=psi;
	MatParams[3]=Kc;
	MatParams[4]=Kphi;
	MatParams[5]=T;
	
	plastic=false;
	inaccurate=0;
	
	// Hardening variables
	aTrial=0.;
	aConvg=0.;

	fSurfaces.resize(2);
	fSurfaces[0]=new DruckerPragerYS(c,phi,Kc,Kphi);
	fSurfaces[1]=new TensionCutOffYS(T);
	
	gSurfaces.resize(2);
	gSurfaces[0]=new DruckerPragerYS(c,psi,Kc,Kphi);
	gSurfaces[1]=new TensionCutOffYS(T);

	// Material tag
	myTag=TAG_MATERIAL_DRUCKER_PRAGER;
}
DruckerPragerNew3::~DruckerPragerNew3()
{
	delete myElastic;
}
MultiaxialMaterial* DruckerPragerNew3::getClone()
{
	// Material parameters
	int myID    = this->getID();
	int elID    = myElastic->getID();
	double c    = MatParams[ 0];
	double phi  = MatParams[ 1];
	double psi  = MatParams[ 2];
	double Kc   = MatParams[ 3];
	double Kphi = MatParams[ 4];
	double T    = MatParams[ 5];
	// Create clone and return
	DruckerPragerNew3* newClone=new DruckerPragerNew3(myID,elID,c,phi,psi,Kc,Kphi,T);
	return newClone;
}


/**
 * Update stresses given a total strain increment.
 * @param De Vector containing total strain increment.
 */ 
void DruckerPragerNew3::setStrain(const Vector& De)
{
	// material properties
	double E = myElastic->getParam(0);
	double nu= myElastic->getParam(1);
	double c    = MatParams[ 0];
	double phi  = MatParams[ 1];
	double psi  = MatParams[ 2];
	double Kc   = MatParams[ 3];
	double Kphi = MatParams[ 4];
	double T    = MatParams[ 5];
	
	// elasticity matrix
	C3(0,0)=  1/E;	C3(0,1)=-nu/E;	C3(0,2)=-nu/E;
	C3(1,0)=-nu/E;	C3(1,1)=  1/E;	C3(1,2)=-nu/E;
	C3(2,0)=-nu/E;	C3(2,1)=-nu/E;	C3(2,2)=  1/E;

	// spectral decomposition
	Vector s(3),De3(3);
	Matrix sV(3,3),eV(3,3);
	aTrial=aConvg;
	eTrial=eTotal+De;
	sTrial=sConvg+(this->getC())*De;
	spectralDecomposition(sTrial,s,sV);
	Vector eTrial3=C3*s;

	Vector snn=sTrial;
	double ann=aTrial;

	//report(eTrial3,"De",8,3);
	//report(s,"De",8,3);

	// Find active surfaces
	//vector<int> activeS;
	//for(unsigned i=0;i<2;i++)
		//if(fSurfaces[i]->getf(s,aTrial)>1e-9)
			//activeS.push_back(i);
	
	vector<bool> activeS;
	for(unsigned i=0;i<2;i++)
	{
        if(fSurfaces[i]->getf(s,aTrial)>1e-9)
		{
            activeS.push_back(true);
		}
		else
		{
			activeS.push_back(false);
		}
	}
	// Elastic (quick return)
	if((activeS[0]==false)&&(activeS[1]==false))
		return;

	// Plastic - start iterations ----------------------------------
	int iter=0;
	Matrix A;
	Vector R;
	Vector x;
	Vector DL(2,0.);
	
	int num_of_iters=50;
	while(iter<num_of_iters)
	{
		++iter;

		// size of the local system
		int nA=0;

		if(activeS[0])		  nA++;		// just Drucker-Prager
		if(activeS[0] && Kc>0) nA++;	// hardening (if Drucker-Prager)
		if(activeS[1])		  nA++;		// tension cut-off
		A.resize(3+nA,3+nA,0.);
		R.resize(3+nA,0.);
		x.resize(3+nA);
		
		// hardening (ok to compute)
		const Vector &v=gSurfaces[0]->getdfds(s,aTrial);
		double dhdq=sqrt(2./3.*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
		double D=Kc*6*cos(phi)/(sqrt(3.)*(3-sin(phi)));

		// build the system
		// [x x x - - -] [x]
		// [x x x - - -] [x]
		// [x x x - - -] [x]
		// [- - - - - -] [-]
		// [- - - - - -] [-]
		// [- - - - - -] [-]
		A.append(C3,0,0,1.0,1.0);
		R.append(C3*s-eTrial3,0,1.0,0.0);
		for(int i=0;i<2;i++) 
		{
			if(activeS[i])
			{
				A.append(DL[i]*(gSurfaces[i]->getd2fdsds(s,aTrial)), 0,   0, 1., 1.);
				R.append(DL[i]*(gSurfaces[i]->getdfds(s,aTrial)),         0, 1., 1.);
			}
		}
		
		// [- - - 0 x x] [-]
		// [- - - 0 x x] [-]
		// [- - - 0 x x] [-]
		// [- - - - - -] [-]
		// [- - - - - -] [-]
		// [- - - - - -] [-]
		int pos=2;
		if(activeS[0] && Kc>0.)
		{
			pos++;
		}
		for(int i=0;i<2;i++) 
		{
			if(activeS[i])
			{
				pos++;
				A.appendCol(gSurfaces[i]->getdfds(s,aTrial),0,pos,1.,0.);
			}
		}

		// [- - - - - -] [-]
		// [- - - - - -] [-]
		// [- - - - - -] [-]
		// [0 0 0 x x 0] [x]
		// [- - - - - -] [-]
		// [- - - - - -] [-]
		if(activeS[0] && Kc>0.)
		{
			A(3,3)=1./D;
			A(3,4)=dhdq;
			R[3]=-aTrial+aConvg+DL[0]*dhdq;
		}

		// [- - - - - -] [-]
		// [- - - - - -] [-]
		// [- - - - - -] [-]
		// [- - - - - -] [-]
		// [x x x 1 0 0] [-]
		// [x x x 0 0 0] [-]
		pos=2;
		if(activeS[0] && Kc>0.)
		{
			A(4,3)=1.;
			pos=3;
		}
		//report(A,"A");	
		for(int i=0;i<2;i++) 
		{
			if(activeS[i])
			{
				pos++;
				A.appendRow(fSurfaces[i]->getdfds(s,aTrial),pos,0,1.,0.);
				R[pos]=fSurfaces[i]->getf(s,aTrial);
			}
		}

		//report(A,"A");			
		//report(R,"R");
		A.solve(x,-R);
		//report(x,"x");

		// update
		s[0]+=x[0];
		s[1]+=x[1];
		s[2]+=x[2];
		pos=2;
		if(activeS[0] && Kc>0.)
		{
			aTrial-=1./D*x[3];
			pos=3;
		}
		for(int i=0;i<2;i++) 
		{
			if(activeS[i])
			{
				pos++;
				DL[i]+=x[pos];
			}
		}

		//cout<<DLambda[0];
		//cout<<R.twonorm()<<endl;
		if(R.twonorm()<1.e-12 && fSurfaces[0]->getf(s,aTrial)<1e-9) break;
	}

	//cout<<aTrial<<endl;
	cout<<iter<<endl;
	if(iter==num_of_iters) report(iter,"error");


	// coordinate transformation
	sTrial[0]=s[0]*sV(0,0)*sV(0,0)+s[1]*sV(1,0)*sV(1,0)+s[2]*sV(2,0)*sV(2,0);
	sTrial[1]=s[0]*sV(0,1)*sV(0,1)+s[1]*sV(1,1)*sV(1,1)+s[2]*sV(2,1)*sV(2,1);
	sTrial[2]=s[0]*sV(0,2)*sV(0,2)+s[1]*sV(1,2)*sV(1,2)+s[2]*sV(2,2)*sV(2,2);
	sTrial[3]=s[0]*sV(0,0)*sV(0,1)+s[1]*sV(1,0)*sV(1,1)+s[2]*sV(2,0)*sV(2,1);
	sTrial[4]=s[0]*sV(0,1)*sV(0,2)+s[1]*sV(1,1)*sV(1,2)+s[2]*sV(2,1)*sV(2,2);
	sTrial[5]=s[0]*sV(0,0)*sV(0,2)+s[1]*sV(1,0)*sV(1,2)+s[2]*sV(2,0)*sV(2,2);

	double Dt =10000.;
	double eta=1000.;
	sTrial=(snn+(Dt/eta)*sTrial)/(1+Dt/eta);
	aTrial=(ann+(Dt/eta)*aTrial)/(1+Dt/eta);

}
/**
 * Commit material state.
 */
void DruckerPragerNew3::commit()
{
	//report(inaccurate);
	inaccurate=0;
	eTotal=eTrial; ///@todo
	sConvg=sTrial;
	aConvg=aTrial;
	this->track();
}
/**
 * Get tangent material matrix.
 * @todo fill it
 * @return A reference to the tangent material matrix.
 */
const Matrix& DruckerPragerNew3::getC()
{
	return myElastic->getC();
}
bool DruckerPragerNew3::isPlastic()
{
	return plastic;
}
/**
 * Add a record to the tracker.
 * If \a myTracker pointer is null (no tracker is added) just return.
 * Otherwise gather info and send them to the tracker.
 * The domain should be already updated!
 */
void DruckerPragerNew3::track()
{
	if(myTracker==0) return;
	ostringstream s;
	s<<"DATA "	<<' ';
//	s<<"sigm "	<<' '<<sConvg;
//	s<<"epst "	<<' '<<eTotal;
//	s<<"epsp "	<<' '<<ePConvg;
//	s<<"epsv "	<<1020<<' '<<eTotal[0]+eTotal[1]+eTotal[2]<<' ';
//	s<<"p "	    <<1020<<' '<<sConvg.p()<<' ';
//	s<<"q "	    <<1020<<' '<<sConvg.q()<<' ';
//	s<<"END "<<' ';
	myTracker->track(pD->getLambda(),pD->getTimeCurr(),s.str());
}
