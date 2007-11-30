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

#include <NemesisDebug.h>

void report(const double d,const char* name,int total,int decimal)
{
	cout<<name<<" =";
	num::print_d(d,total,decimal);
	cout<<endl;
}
void report(const Matrix& m,const char* name,int total,int decimal)
{
	cout<<name<<" [ "<<m.rows()<<","<<m.cols()<<" ] ="<<endl;
	for(int i=0;i<m.rows();i++)
	{
		for(int j=0;j<m.cols();j++)
			num::print_d(m(i,j),total,decimal);
		cout<<endl;
	}
}
void report(const Vector& v,const char* name,int total,int decimal)
{
	cout<<name<<" [ "<<v.size()<<" ] =";
	for(int i=0;i<v.size();i++)
		num::print_d(v[i],total,decimal);
	cout<<endl;
}


void add(Matrix& K,int row,int col,const Matrix& B1,const Matrix& C,const Matrix B2,double c1,double c0)
{
	int m=B1.rows();
	int n=B2.cols();
	int pos=row*K.cols()+col;
	int colK=K.cols();
	double* pK =K.data();
	double* pB1=B1.data();
	double* pB2=B2.data();
	double* pC=C.data();
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			pK[pos+i*colK+j]*=c0;
	for(int i=0;i<n;i++)
		for(int k=0;k<m;k++)
			for(int l=0;l<m;l++)
				for(int j=0;j<n;j++)
					pK[pos+i*colK+j]+=c1*pB1[k*n+i]*pC[k*m+l]*pB2[l*n+j];
}
void add(Vector& R,int row,const Matrix& BT,const Vector& V,double c1,double c0)
{
	int m=BT.rows();
	int n=BT.cols();
	double* pV =V.data();
	double* pR =R.data();
	double* pBT=BT.data();
	for(int i=0;i<n;i++)
		pR[i]*=c0;
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			pR[row+i]+=c1*pBT[j*n+i]*pV[j];
}
void add2(Vector& R,int row,const Matrix& BT,const Vector& V,double c1,double c0)
{
	int m=BT.rows();
	int n=BT.cols();
	double* pV =V.data();
	double* pR =R.data();
	double* pBT=BT.data();
	for(int i=0;i<m;i++)
		pR[i]*=c0;
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			pR[i]+=c1*pBT[i*n+j]*pV[row+j];
}
//*****************************************************************************
//
//*****************************************************************************
void add_BTCB(Matrix& K,int row,int col,const int* perm,const Matrix& B1,const Matrix& C,const Matrix B2,double c1,double c0)
{
	int m=B1.rows();
	int n=B2.cols();
	int pos=row*K.cols()+col;
	int colK=K.cols();
	int colC=C.cols();
	double* pK =K.data();
	double* pB1=B1.data();
	double* pB2=B2.data();
	double* pC=C.data();
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			pK[pos+i*colK+j]*=c0;
	for(int i=0;i<n;i++)
		for(int k=0;k<m;k++)
			for(int l=0;l<m;l++)
				for(int j=0;j<n;j++)
					pK[pos+i*colK+j]+=c1*pB1[k*n+i]*pC[perm[k]*colC+perm[l]]*pB2[l*n+j];
}
void add_BTv (Vector& R,int row,const int* perm,const Matrix& B,const Vector& v,double c1,double c0)
{
	int m=B.rows();
	int n=B.cols();
	double* pV=v.data();
	double* pR=R.data();
	double* pB=B.data();
	for(int i=0;i<n;i++)
		pR[i]*=c0;
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			pR[row+i]+=c1*pB[j*n+i]*pV[perm[j]];
}
void add_Bv  (Vector& R,int row,const int* perm,const Matrix& B,const Vector& v,double c1,double c0)
{
	int m=B.rows();
	int n=B.cols();
	double* pV=v.data();
	double* pR=R.data();
	double* pB=B.data();
	for(int i=0;i<m;i++)
		pR[i]*=c0;
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			pR[perm[i]]+=c1*pB[i*n+j]*pV[row+j];
}