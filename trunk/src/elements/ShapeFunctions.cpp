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

// Included files
#include <ShapeFunctions.h>

void shape4(const Matrix& x,double shp[4][3][4],double detJ[4])
{
	const double dsq3=0.577350269189626;
	static const double gCrds[4][2]={	{-dsq3,-dsq3},
										{+dsq3,-dsq3},
									    {+dsq3,+dsq3},
										{-dsq3,+dsq3}};
	for(int k=0;k<4;k++)
	{
		double m1=1-gCrds[k][0];		// 1-xi
		double p1=1+gCrds[k][0];		// 1+xi
		double m2=1-gCrds[k][1];		// 1-eta
		double p2=1+gCrds[k][1];		// 1+eta
		
		shp[0][0][k]= 0.25*m1*m2;		// N1
		shp[1][0][k]= 0.25*p1*m2;		// N2
		shp[2][0][k]= 0.25*p1*p2;		// N3
		shp[3][0][k]= 0.25*m1*p2;		// N4

		shp[0][1][k]=-0.25*m2;			// N1,xi
		shp[1][1][k]=+0.25*m2;			// N2,xi
		shp[2][1][k]=+0.25*p2;			// N3,xi
		shp[3][1][k]=-0.25*p2;			// N4,xi

		shp[0][2][k]=-0.25*m1;			// N1,eta
		shp[1][2][k]=-0.25*p1;			// N2,eta
		shp[2][2][k]=+0.25*p1;			// N3,eta
		shp[3][2][k]=+0.25*m1;			// N4,eta
	
		// Jacobian
		static double J[2][2];
		for(int i=0;i<2;i++)
		{
			J[0][i]=0.;
			J[1][i]=0.;
			for(int j=0;j<4;j++)
			{
				J[0][i]+=shp[j][1][k]*x(j,i);
				J[1][i]+=shp[j][2][k]*x(j,i);
			}
		}

		// Determinant J
        detJ[k]=J[0][0]*J[1][1]-J[0][1]*J[1][0];
        double ddetJ=1./detJ[k];
		
		// Inverse J
		double d;
		d=J[0][0];
		J[0][0]= ddetJ*J[1][1];
		J[1][1]= ddetJ*d;
		d=J[1][0];
		J[1][0]=-ddetJ*J[0][1];
		J[0][1]=-ddetJ*d;

		// dNdx
		for(int i=0;i<4;i++)
		{
			double d1=shp[i][1][k]*J[0][0]+shp[i][2][k]*J[0][1];
			double d2=shp[i][1][k]*J[1][0]+shp[i][2][k]*J[1][1];
			shp[i][1][k]=d1;
			shp[i][2][k]=d2;
		}
	}
}

void shape8(const Matrix& x,double shp[8][4][8],double detJ[8])
{
	const double dsq3=0.577350269189626;
	static const double gCrds[8][3]={	{-dsq3,-dsq3,-dsq3},
										{+dsq3,-dsq3,-dsq3},
										{+dsq3,+dsq3,-dsq3},
										{-dsq3,+dsq3,-dsq3},
										{-dsq3,-dsq3,+dsq3},
										{+dsq3,-dsq3,+dsq3},
										{+dsq3,+dsq3,+dsq3},
										{-dsq3,+dsq3,+dsq3}};
	for(int k=0;k<8;k++)
	{
		double m1=1-gCrds[k][0];		// 1-xi
		double p1=1+gCrds[k][0];		// 1+xi
		double m2=1-gCrds[k][1];		// 1-eta
		double p2=1+gCrds[k][1];		// 1+eta
		double m3=1-gCrds[k][2];		// 1-zeta
		double p3=1+gCrds[k][2];		// 1+zeta
		
		shp[0][0][k]= 0.125*m1*m2*m3;	// N1
		shp[1][0][k]= 0.125*p1*m2*m3;	// N2
		shp[2][0][k]= 0.125*p1*p2*m3;	// N3
		shp[3][0][k]= 0.125*m1*p2*m3;	// N4
		shp[4][0][k]= 0.125*m1*m2*p3;	// N5
		shp[5][0][k]= 0.125*p1*m2*p3;	// N6
		shp[6][0][k]= 0.125*p1*p2*p3;	// N7
		shp[7][0][k]= 0.125*m1*p2*p3;	// N8

		shp[0][1][k]=-0.125*m2*m3;		// N1,xi
		shp[1][1][k]=+0.125*m2*m3;		// N2,xi
		shp[2][1][k]=+0.125*p2*m3;		// N3,xi
		shp[3][1][k]=-0.125*p2*m3;		// N4,xi
		shp[4][1][k]=-0.125*m2*p3;		// N5,xi
		shp[5][1][k]=+0.125*m2*p3;		// N6,xi
		shp[6][1][k]=+0.125*p2*p3;		// N7,xi
		shp[7][1][k]=-0.125*p2*p3;		// N8,xi

		shp[0][2][k]=-0.125*m1*m3;		// N1,eta
		shp[1][2][k]=-0.125*p1*m3;		// N2,eta
		shp[2][2][k]=+0.125*p1*m3;		// N3,eta
		shp[3][2][k]=+0.125*m1*m3;		// N4,eta
		shp[4][2][k]=-0.125*m1*p3;		// N5,eta
		shp[5][2][k]=-0.125*p1*p3;		// N6,eta
		shp[6][2][k]=+0.125*p1*p3;		// N7,eta
		shp[7][2][k]=+0.125*m1*p3;		// N8,eta

		shp[0][3][k]=-0.125*m1*m2;		// N1,zeta
		shp[1][3][k]=-0.125*p1*m2;		// N2,zeta
		shp[2][3][k]=-0.125*p1*p2;		// N3,zeta
		shp[3][3][k]=-0.125*m1*p2;		// N4,zeta
		shp[4][3][k]=+0.125*m1*m2;		// N5,zeta
		shp[5][3][k]=+0.125*p1*m2;		// N6,zeta
		shp[6][3][k]=+0.125*p1*p2;		// N7,zeta
		shp[7][3][k]=+0.125*m1*p2;		// N8,zeta
	
		// Jacobian
		static double J[3][3];
		for(int i=0;i<3;i++)
		{
			J[0][i]=0.;
			J[1][i]=0.;
			J[2][i]=0.;
			for(int j=0;j<8;j++)
			{
				J[0][i]+=shp[j][1][k]*x(j,i);
				J[1][i]+=shp[j][2][k]*x(j,i);
				J[2][i]+=shp[j][3][k]*x(j,i);
			}
		}

		static double cof[3][3];
		cof[0][0]=J[1][1]*J[2][2]-J[2][1]*J[1][2];
		cof[0][1]=J[2][1]*J[0][2]-J[0][1]*J[2][2];
		cof[0][2]=J[0][1]*J[1][2]-J[1][1]*J[0][2];

		cof[1][0]=J[1][2]*J[2][0]-J[2][2]*J[1][0];
		cof[1][1]=J[2][2]*J[0][0]-J[0][2]*J[2][0];
		cof[1][2]=J[0][2]*J[1][0]-J[1][2]*J[0][0];

		cof[2][0]=J[1][0]*J[2][1]-J[2][0]*J[1][1];
		cof[2][1]=J[2][0]*J[0][1]-J[0][0]*J[2][1];
		cof[2][2]=J[0][0]*J[1][1]-J[1][0]*J[0][1];

		// Determinant J
        detJ[k]=J[0][0]*cof[0][0]+J[0][1]*cof[1][0]+J[2][0]*cof[0][2];
		
		// Inverse J
		double ddetJ=1./detJ[k];
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				J[i][j]=cof[i][j]*ddetJ;
		
		// dNdx
		for(int i=0;i<8;i++)
		{
			double d1=shp[i][1][k]*J[0][0]+shp[i][2][k]*J[0][1]+shp[i][3][k]*J[0][2];
			double d2=shp[i][1][k]*J[1][0]+shp[i][2][k]*J[1][1]+shp[i][3][k]*J[1][2];
			double d3=shp[i][1][k]*J[2][0]+shp[i][2][k]*J[2][1]+shp[i][3][k]*J[2][2];
			shp[i][1][k]=d1;
			shp[i][2][k]=d2;
			shp[i][3][k]=d3;
		}
	}
}


