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

#include "main/nemesis_debug.h"
#include "material/tresca.h"

Matrix Tresca::C(6,6,0.);
Matrix Tresca::C3(3,3,0.);

Tresca::Tresca()
{
}
Tresca::Tresca(int ID,int elasticID,double cu,double kx,double ky,double kz)
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
	MatParams[0]=cu;
	MatParams[1]=kx;
	MatParams[2]=ky;
	MatParams[3]=kz;
	
	// Material state
//	ePTrial.resize(6,0.);	ePConvg.resize(6,0.);
//	qTrial.resize(6,0.);	qConvg.resize(6,0.);
//	aTrial=0.;				aConvg=0.;
	
	plastic=false;
	inaccurate=0;
}
Tresca::~Tresca()
{
	delete myElastic;
}
MultiaxialMaterial* Tresca::getClone()
{
	// Material parameters
	int myID    = this->getID();
	int elID    = myElastic->getID();
	double cu   = MatParams[ 0];
	double kx   = MatParams[ 1];
	double ky   = MatParams[ 2];
	double kz   = MatParams[ 3];
	// Create clone and return
	Tresca* newClone=new Tresca(myID,elID,cu,kx,ky,kz);
	return newClone;
}
/**
 * Update stresses given a total strain increment.
 * @param De Vector containing total strain increment.
 */ 
void Tresca::setStrain(const Vector& De)
{
	std::vector<Vector> df(3);
	df[0].resize(3);	df[1].resize(3);	df[2].resize(3);
	df[0][0]= 1.0;		df[0][1]= 0.0;		df[0][2]=-1.0;	
	df[1][0]= 1.0;		df[1][1]=-1.0;		df[1][2]= 0.0;	
	df[2][0]= 0.0;		df[2][1]= 1.0;		df[2][2]=-1.0;	
	
	double E = myElastic->getParam(0);
	double nu= myElastic->getParam(1);
	double kx   = MatParams[ 1];
	double ky   = MatParams[ 2];
	double kz   = MatParams[ 3];
	double cu= MatParams[ 0]+(kx*x+ky*y+kz*z);
	C3(0,0)=  1/E;	C3(0,1)=-nu/E;	C3(0,2)=-nu/E;
	C3(1,0)=-nu/E;	C3(1,1)=  1/E;	C3(1,2)=-nu/E;
	C3(2,0)=-nu/E;	C3(2,1)=-nu/E;	C3(2,2)=  1/E;

	static Vector s(3);
	static Matrix sV(3,3);
	eTrial=eTotal+De;
	sTrial=sConvg+(this->getC())*De;
	spectralDecomposition(sTrial,s,sV);
	//Vector tempS(3);
	//tempS=sTrial;

	static Vector f(3);
	f[0]=s[0]-s[2]-2.*cu;
	f[1]=s[0]-s[1]-2.*cu;
	f[2]=s[1]-s[2]-2.*cu;
	//report(f,"Yield fun",true,12);
	double theta=sTrial.theta();
	//report(2.*sqrt(sTrial.J2())*cos(theta)-cu,"f");

	std::vector<int> active;
	for(unsigned i=0;i<3;i++)
		if(f[i]>0.) active.push_back(i);

	// Elastic case
	if(active.size()==0) return;
	plastic=true;

	// Plastic case
	static Matrix A;
	static Vector x;
	static Vector R;

		active.clear();
		for(unsigned i=0;i<3;i++)
			if(f[i]>0.) active.push_back(i);
	for(int k=0;k<4;k++)
	{
		A.resize(3+active.size(),3+active.size(),0.);
		x.resize(3+active.size());
		R.resize(3+active.size());
		R.clear();
		
		A.append(C3,0,0);
		for(unsigned i=0;i<active.size();i++)
		{
			A.appendCol(df[active[i]],  0,3+i);
			A.appendRow(df[active[i]],3+i,  0);
			R[3+i]=-f[active[i]];
		}
		
		// solve
		//report(A,"A");
		//report(R);
		A.solve(x,R);
		//report(x);

		// check
		bool restart=false;
		for(unsigned i=0;i<active.size();i++)
		{
			if(x[3+i]<0.)
			{
				active.erase(active.begin()+i,active.begin()+i+1);
				restart=true;
			}
		}
		if(restart) continue;
		
		// update
		for(int i=0;i<3;i++) s[i]+=x[i];
		break;
	}

//	cout<<s[0]-s[2]-cu<<endl;
//	cout<<s[0]-s[1]-cu<<endl;
//	cout<<s[1]-s[2]-cu<<endl;
	// coordinate transformation
	sTrial[0]=s[0]*sV(0,0)*sV(0,0)+s[1]*sV(1,0)*sV(1,0)+s[2]*sV(2,0)*sV(2,0);
	sTrial[1]=s[0]*sV(0,1)*sV(0,1)+s[1]*sV(1,1)*sV(1,1)+s[2]*sV(2,1)*sV(2,1);
	sTrial[2]=s[0]*sV(0,2)*sV(0,2)+s[1]*sV(1,2)*sV(1,2)+s[2]*sV(2,2)*sV(2,2);
	sTrial[3]=s[0]*sV(0,0)*sV(0,1)+s[1]*sV(1,0)*sV(1,1)+s[2]*sV(2,0)*sV(2,1);
	sTrial[4]=s[0]*sV(0,1)*sV(0,2)+s[1]*sV(1,1)*sV(1,2)+s[2]*sV(2,1)*sV(2,2);
	sTrial[5]=s[0]*sV(0,0)*sV(0,2)+s[1]*sV(1,0)*sV(1,2)+s[2]*sV(2,0)*sV(2,2);

	theta=sTrial.theta();
	if((2.*sqrt(sTrial.J2())*cos(theta)-2*cu)>1e-6) 
	{
		//inaccurate++;
		//cout<<"error : "<<2.*sqrt(sTrial.J2())*cos(theta)-2*cu<<endl;
		//report(tempS,"sTrial");
		//report(f,    "yield");
		//report(x,    "solns",20,15);
        //exit(0);
	}
	//for(unsigned i=0;i<active.size();i++)
	//	cout<<active[i]<<endl;
}
/**
 * Commit material state.
 */
void Tresca::commit()
{
	//report(inaccurate);
	inaccurate=0;
	eTotal=eTrial; ///@todo
	sConvg=sTrial;
	this->track();
}
/**
 * Get tangent material matrix.
 * @todo fill it
 * @return A reference to the tangent material matrix.
 */
const Matrix& Tresca::getC()
{
	return myElastic->getC();
}
bool Tresca::isPlastic()
{
	return plastic;
}
/**
 * Add a record to the tracker.
 * If \a myTracker pointer is null (no tracker is added) just return.
 * Otherwise gather info and send them to the tracker.
 * The domain should be already updated!
 */
void Tresca::track()
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
