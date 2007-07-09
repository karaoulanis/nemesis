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

#include <MultiaxialElastoPlastic.h>
#include <NemesisDebug.h>

Matrix MultiaxialElastoPlastic::C(6,6,0.);

MultiaxialElastoPlastic::MultiaxialElastoPlastic()
{
}
MultiaxialElastoPlastic::MultiaxialElastoPlastic(int ID,int elasticID)
:MultiaxialMaterial(ID,0.,0.)
{
	// Get the elastic part
	Material* p=pD->get<Material>(pD->getMaterials(),elasticID);
	if(p->getTag()!=TAG_MATERIAL_MULTIAXIAL_ELASTIC)
		throw SException("[nemesis:%d] %s",9999,"Multiaxial elastic material expected.");
	myElastic=static_cast<MultiaxialMaterial*>(p)->getClone();
	MatParams[30]=myElastic->getParam(30);
	MatParams[31]=myElastic->getParam(31);

	// Material state
	eTotal.resize(6,0.);
	ePTrial.resize(6,0.);	ePConvg.resize(6,0.);
	qTrial.resize(6,0.);	qConvg.resize(6,0.);
	aTrial.resize(6,0.);	aConvg.resize(6,0.);
	plastic=false;
}
MultiaxialElastoPlastic::~MultiaxialElastoPlastic()
{
	delete myElastic;
	Containers::vector_delete(fSurfaces);
	Containers::vector_delete(gSurfaces);
}
/**
 * Update stresses given a total strain increment.
 * @param De Vector containing total strain increment.
 */ 
void MultiaxialElastoPlastic::setStrain(const Vector& De)
{
	//this->returnMapSYS(De);
	//this->returnMapTest(De);
	if(fSurfaces.size()==1)	this->returnMapSYS(De);
	//if(fSurfaces.size()==1)	this->returnMapTest(De);
	else					this->returnMapMYS(De);
}
/**
 * Commit material state.
 */
void MultiaxialElastoPlastic::commit()
{
	ePConvg=ePTrial;
	sConvg=sTrial;
	qConvg=qTrial;
	aConvg=aTrial;
}
void MultiaxialElastoPlastic::returnMapTest(const Vector& De)
{
	static LogFile log("TYS.log");

	int nIter=20;
	double tol1=1e-6;
	double tol2=1e-6;
	Vector R(7,0.);
	Matrix Cel=myElastic->getC();
	Matrix invCel=Inverse(Cel);
	Matrix A(7,7,0.);
	ePTrial=ePConvg;
	aTrial[0]=aConvg[0];
	double dg=0;
	Surface* fS=fSurfaces[0];
	double eta=1000.;
	double dt=pD->getTimeIncr();
	
	//=========================================================================
	// Step 1: Compute trial stress
	//=========================================================================
	sTrial=sConvg+Cel*De;
	
	//=========================================================================
	// Step 2: Check for plastic process
	//=========================================================================
	if(fS->get_f(sTrial)<=0) return;
	Counters::c2=0;
 	for(int k=0;k<nIter;k++)
	{
		// Residual Vector 
		Vector tmp(6,0.);
		tmp=De;
		tmp-=invCel*(sTrial-sConvg);
		//if(k==0) log.write(tmp);
		tmp-=dg*(fS->get_dfds(sTrial));
		R.clear();
		R.append(tmp,0);
		//R[6]=-(fS->get_f(sTrial))+(eta/dt)*dg;
		R[6]=-(fS->get_f(sTrial));
		//log<< ++Counters::c2<<'\n';
		//log.write(sConvg);
		//log.write(sTrial);
		//log.write(fS->get_dfds(sTrial));
		//log.write(invCel*(sTrial-sConvg));
		// Convergence check
		if(tmp.twonorm()<tol1 && fS->get_f(sTrial)<tol2)
		{
			aTrial[0]+=dg;
			break;
		}
		// Jacobian
		A.clear();
		A.append(invCel+dg*(fS->get_df2dss(sTrial)),0,0);
		for(int i=0;i<6;i++) A(6,i)=(fS->get_dfds(sTrial))[i];
		for(int i=0;i<6;i++) A(i,6)=(fS->get_dfds(sTrial))[i];
		//A(6,6)=-eta/dt;
		//for(int i=0;i<7;i++)
		//{
		//	for(int j=0;j<7;j++) cout<<A(i,j)<<"\t\t";
		//	cout<<endl;
		//}
		//cout<<endl;
		// Solution and update
		
		//log.write(A);
		//if(++Counters::c1==60) log.write(A);
		//for(int i=0;i<6;i++) for(int j=0;j<6;j++) log<<A(i,j)<<" "; log<<'\n'; 
		//log.write(sTrial); 

		Vector x(7,0.);
		A.solve(x,R);
		for(int i=0;i<6;i++) sTrial[i]+=x[i];
		dg+=x[6];
		log<<dg<<'\n';
	}
}
/**
 * Single surface return mapping.
 * @param De Vector containing total strain increment.
 */
void MultiaxialElastoPlastic::returnMapSYS(const Vector& De)
{
	static LogFile log("SYS.log");
	//=========================================================================
	// Setup
	//=========================================================================
	int nIter=25;
	double tol1=1e-6;
	double tol2=1e-6;
	double ddg;
	Vector R(6,0.);
	Vector dEp(6,0.);
	Matrix Cel   =myElastic->getC();
	Matrix invCel=Inverse(Cel);
	
	ePTrial=ePConvg;
	aTrial[0]=0.;
	//aTrial[0]=aConvg[0];

	//double& dg=aTrial[0];
	double dg=0.;
	Surface* fS=fSurfaces[0];
	Surface* gS=gSurfaces[0];
	
	Vector enn=invCel*sConvg+ePConvg+De;

	//=========================================================================
	// Step 1: Compute trial stress
	//=========================================================================
	sTrial=sConvg+Cel*De;
	Vector ss=sTrial;
 
	//=========================================================================
	// Step 2: Check for plastic process
	//=========================================================================
	if(fS->get_f(sTrial)<=0) return;

	Counters::c2=0;
	for(int k=0;k<nIter;k++)
	{
		//=====================================================================
		// Step 3: Evaluate flow rule, hardening law and residuals
		//=====================================================================
		sTrial=Cel*(enn-ePTrial);
		R=-ePTrial+ePConvg+dg*(fS->get_dfds(sTrial));
		//if(k==0) log.write(-ePTrial+ePConvg);
		//log<< ++Counters::c2<<'\n';
		//log.write(sConvg);
		//log.write(sTrial);
		//log.write(fS->get_dfds(sTrial));
		//log.write(-ePTrial+ePConvg);
		//for(int i=0;i<6;i++) log<<R[i]<<" "; log<<'\n';

		//=====================================================================
		// Step 4: Check convergence
		//=====================================================================
		if(R.twonorm()<tol1 && fS->get_f(sTrial)<tol2) break;

		//=====================================================================
		// Step 5: Compute elastic moduli and consistent tangent moduli
		//=====================================================================
		// Matrix A 
		Matrix A(6,6,0.);
		A.append(invCel+dg*(fS->get_df2dss(sTrial)),0,0);
		//log.write(A);
		//if(++Counters::c1==60) log.write(A);
		//for(int i=0;i<6;i++) for(int j=0;j<6;j++) log<<A(i,j)<<" "; log<<'\n'; 
		//for(int i=0;i<6;i++) log.write(sTrial); 
		//log.write(sTrial); 
		A=Inverse(A);

		//=====================================================================
		// Step 6: Obtain increment to consistency parameter
		//=====================================================================
		ddg=(fS->get_f(sTrial)-fS->get_dfds(sTrial)*(A*R))/(fS->get_dfds(sTrial)*(A*(fS->get_dfds(sTrial))));

		//=====================================================================
		// Step 7: Obtain incremental plastic strains and internal variables
		//=====================================================================
		dEp=invCel*A*(R+ddg*(fS->get_dfds(sTrial)));

		//=====================================================================
		// Step 8: Update
		//=====================================================================
		ePTrial+=dEp;
		dg+=ddg;
		log<<dg<<'\n';
		//cout<<k<<" "<<dg<<" "<<R.twonorm()<<endl;
	}
	//cout<<endl;
	double dt=pD->getTimeIncr();
	double eta=1000.;
	//sTrial=((sConvg+Cel*De)+(dt/eta)*(sConvg+Cel*De))/(1+(dt/eta));
	//sTrial=(ss+(dt/eta)*sTrial)/(1+(dt/eta));
}
/**
 * Multiple surface return mapping.
 * @param De Vector containing total strain increment.
 */
void MultiaxialElastoPlastic::returnMapMYS(const Vector& De)
{
	//=========================================================================
	// Setup
	//=========================================================================
	int nIter=30;
	double tol1=1e-6;
	double tol2=1e-6;
	Vector ddg(12,0.);
	Vector R(6,0.);
	Vector dEp(6,0.);
	Matrix Cel   =myElastic->getC();
	Matrix invCel=Inverse(Cel);
	
	ePTrial=ePConvg;//todo:CHECK!!!!
	aTrial=aConvg;	//todo:CHECK!!!!
	aTrial.clear();	//todo:CHECK!!!!
	Vector& dg=aTrial;
	Vector enn=invCel*sConvg+ePConvg+De;
	
	//=========================================================================
	// Step 1: Compute trial stress
	//=========================================================================
	sTrial=sConvg+Cel*De;
	Vector ss=sTrial;
	plastic=false;
   
	//=========================================================================
	// Step 2: Check for plastic process
	//=========================================================================
	int nActiveSurfaces=0;
	for(unsigned i=0;i<fSurfaces.size();i++)
		if(fSurfaces[i]->get_f(sTrial)>0)
		{
			fSurfaces[i]->setActive(true);
			nActiveSurfaces++;
		}
	if(nActiveSurfaces==0) return;
	
	plastic=true;
	for(int k=0;k<nIter;k++)
	{
		//=====================================================================
		// Step 3: Evaluate flow rule, hardening law and residuals
		//=====================================================================
		sTrial=Cel*(enn-ePTrial);
		R=-ePTrial+ePConvg;
		for(unsigned b=0;b<fSurfaces.size();b++)
			if(fSurfaces[b]->isActive())
				R+=dg[b]*(fSurfaces[b]->get_dfds(sTrial));
		//cout<<ePTrial<<endl;

		//=====================================================================
		// Step 4: Check convergence
		//=====================================================================
		bool converged=true;
		if(R.twonorm()>tol1) converged=false;
		for(unsigned a=0;a<fSurfaces.size();a++)
			if(fSurfaces[a]->isActive())
				if(abs(fSurfaces[a]->get_f(sTrial))>tol2) converged=false;

		if(converged) break;

		//=====================================================================
		// Step 5: Compute elastic moduli and consistent tangent moduli
		//=====================================================================
		// Matrix A 
		Matrix A(6,6,0.);
		A.append(invCel,0,0);
		for(unsigned b=0;b<fSurfaces.size();b++)
			if(fSurfaces[b]->isActive())
				A+=dg[b]*(fSurfaces[b]->get_df2dss(sTrial));
		A=Inverse(A);

		// Matrix G
		Matrix Gab(4,4,0.);
		Gab(0,0)=1.0;Gab(1,1)=1.0;Gab(2,2)=1.0;Gab(3,3)=1.0;
		for(unsigned a=0;a<fSurfaces.size();a++)
			for(unsigned b=0;b<fSurfaces.size();b++)
				if(fSurfaces[a]->isActive() && fSurfaces[a]->isActive())
				{
					Vector va(6,0.);
					Vector vb(6,0.);
					va=fSurfaces[a]->get_dfds(sTrial);
					vb=fSurfaces[b]->get_dfds(sTrial);
	                Gab(a,b)=(va*(A*vb));
				}
		Gab=Inverse(Gab);

		//=====================================================================
		// Step 6: Obtain increment to consistency parameter
		//=====================================================================
		ddg.clear();
		for(unsigned a=0;a<fSurfaces.size();a++)
			for(unsigned b=0;b<fSurfaces.size();b++)
				if(fSurfaces[a]->isActive() && fSurfaces[b]->isActive())
				{
					Vector vb=fSurfaces[b]->get_dfds(sTrial);
					double fb=fSurfaces[b]->get_f(sTrial);
					ddg[a]+=Gab(a,b)*(fb-vb*(A*R));
				}

		bool reset=false;
		for(unsigned a=0;a<fSurfaces.size();a++)
			if(dg[a]+ddg[a]<0.)
				if(fSurfaces[a]->isActive())
				{
					fSurfaces[a]->setActive(false);
					dg[a]=0.;
					ddg[a]=0.;
					reset=true;
				}
		if(reset) continue;

		//=====================================================================
		// Step 7: Obtain incremental plastic strains and internal variables
		//=====================================================================
		static Vector temp(6);
		temp.clear();
		for(unsigned b=0;b<fSurfaces.size();b++)
			if(fSurfaces[b]->isActive())
				temp+=ddg[b]*(fSurfaces[b]->get_dfds(sTrial));

		//=====================================================================
		// Step 8: Update
		//=====================================================================
		dEp=invCel*A*(R+temp);
		ePTrial+=dEp;
		for(unsigned a=0;a<fSurfaces.size();a++)
			if(fSurfaces[a]->isActive()) dg[a]+=ddg[a];
	}
	if(k==nIter)
	{
//		cout.precision(12);
//		cout<<sConvg+Cel*De<<endl;
//		cout<<De<<endl;
//		cout<<sConvg<<endl;
//		cout<<ePConvg<<endl;
//		cout<<ePTrial<<endl;
//		cout<<"PASSED : "<<nActiveSurfaces<<endl;
	}
	//cout<<endl;
	double dt=pD->getTimeIncr();
	double eta=1000.;
	//sTrial=((sConvg+Cel*De)+(dt/eta)*(sConvg+Cel*De))/(1+(dt/eta));
	//sTrial=(ss+(dt/eta)*sTrial)/(1+(dt/eta));
}
/**
 * Get tangent material matrix.
 * @todo fill it
 * @return A reference to the tangent material matrix.
 */
const Matrix& MultiaxialElastoPlastic::getC()
{
	return myElastic->getC();
}
