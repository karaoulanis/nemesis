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

#include "elements/brick8d.h"
#include "main/NemesisDebug.h"

Brick8d::Brick8d()
{
}
Brick8d::Brick8d(int ID,
					   int Node_1,int Node_2,int Node_3,int Node_4,
					   int Node_5,int Node_6,int Node_7,int Node_8,
					   int matID)
:Brick8(ID,Node_1,Node_2,Node_3,Node_4,Node_5,Node_6,Node_7,Node_8,matID)
{
	myTag=TAG_ELEM_BRICK_8_DISP;
}
Brick8d::~Brick8d()
{
}
void Brick8d::getB(Matrix& B,int node,int gPoint)
{
	// Resize B
	B.resize(6,3);
	double B1=shp[node][1][gPoint];
	double B2=shp[node][2][gPoint];
	double B3=shp[node][3][gPoint];

	B(0,0)=B1;	B(0,1)=0.;	B(0,2)=0.;
	B(1,0)=0.;	B(1,1)=B2;	B(1,2)=0.;
	B(2,0)=0.;	B(2,1)=0.;	B(2,2)=B3;
	B(3,0)=B2;	B(3,1)=B1;	B(3,2)=0.;
	B(4,0)=0.;	B(4,1)=B3;	B(4,2)=B2;
	B(5,0)=B3;	B(5,1)=0.;	B(5,2)=B1;
}
