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

#ifndef _TIMOSHENKO2D_H
#define _TIMOSHENKO2D_H

#include "elements/element.h"
#include "material/uniaxial_material.h"

class Timoshenko2d: public Element
{
protected:
	double cosX[2];
	int mySecID;
	CrossSection* mySection;
	double L;
	UniaxialMaterial* myUniMaterial;
	static const double GaussCoords[4][4];
	static const double GaussWeights[4][4];
	int gPoints;
public:
	Timoshenko2d();
	Timoshenko2d(int ID,int Node_1,int Node_2,int matID,int secID,int rule);	
	Timoshenko2d(int ID,int Node_1,int Node_2,int Node_3,int matID,int secID,int rule);	
	~Timoshenko2d();

	const Matrix& getK();
    const Matrix& getM();
	const Vector& getR();
	const Vector& getRgrad();

	bool checkIfAllows(FEObject* /*f*/)		{return true;}
	void update()						{return;}
	void commit()						{return;}
	void recoverStresses();
	void shapeFunctions(int n,double xi,double &N, double &dN);
};
#endif
