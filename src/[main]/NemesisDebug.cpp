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
void report(const Vector& v,const char* name,bool transpose,int total,int decimal)
{
	cout<<name<<" [ "<<v.size()<<" ] =";
	if(!transpose) cout<<endl;
	for(int i=0;i<v.size();i++)
	{
		num::print_d(v[i],total,decimal);
		if(!transpose) cout<<endl;
	}
	if(transpose) cout<<endl;
}