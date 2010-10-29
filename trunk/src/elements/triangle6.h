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

#ifndef _TRIANGLE6_H
#define _TRIANGLE6_H

#include "elements/element.h"
#include "material/matpoint.h"

class Triangle6: public Element
{
protected:
	double a1,a2,a3;
	double b1,b2,b3;
	double c1,c2,c3;
	double A;
	std::vector<MatPoint*> myMatPoints;
public:
	// Constructors and Destructor
	Triangle6();
	Triangle6(int ID,int Node_1,int Node_2,int Node_3,
		int Node_4,int Node_5,int Node_6,int matID);
	~Triangle6();

	const Matrix& getK();
	const Matrix& getM();
	const Vector& getR();

	void update();
	void commit();
	
	void getShapeFunctions(MatPoint* pMatPoint,Matrix& N,double& detJ);
	void addInitialStresses(InitialStresses* pInitialStresses);
	bool checkIfAllows(FEObject* f);
	void recoverStresses();
	int getnPlasticPoints();
};

#endif
