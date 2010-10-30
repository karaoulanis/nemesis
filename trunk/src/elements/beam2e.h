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

#ifndef NEMESIS_ELEMENTS_BEAM2E_H_
#define NEMESIS_ELEMENTS_BEAM2E_H_

#include "elements/element.h"
#include "material/uniaxial_material.h"

class Beam2e: public Element
{
protected:
	double cosX[2];
	int mySecID;
	CrossSection* mySection;
	double L;
	UniaxialMaterial* myUniMaterial;
public:
	Beam2e();
	Beam2e(int ID,int Node_1,int Node_2,int matID,int secID);	
	~Beam2e();

	const Matrix& getK();
    const Matrix& getM();
	const Vector& getR();
	const Vector& getRgrad();

	bool checkIfAllows(FEObject* /*f*/)		{return true;}
	void update()						{return;}
	void commit()						{return;}
	void recoverStresses();
};
#endif //NEMESIS_ELEMENTS_BEAM2E_H_
