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

#ifndef _BRICK8B_H
#define _BRICK8B_H

#include <Element.h>
#include <MatPoint.h>

class Brick8b: public Element
{
protected:
	static Vector N;
	static Matrix B[8];
	static double detJ;
	static bool first;
	
	Vector Bb1;
	Vector Bb2;
	Vector Bb3;
	std::vector<MatPoint*> myMatPoints;
public:
	// Constructors and Destructor
	Brick8b();
	Brick8b(int ID,
				int Node_1,int Node_2,int Node_3,int Node_4,	
				int Node_5,int Node_6,int Node_7,int Node_8,
				int matID);	
	~Brick8b();

	const Matrix& getK();
    const Matrix& getM();
	const Vector& getR();

	void update();
	void commit();
	
	void shapeFunctions(MatPoint* pMatPoint);
	bool checkIfAllows(FEObject* f);
	void addInitialStresses(InitialStresses* pInitialStresses);	
	void recoverStresses();
	const int getnPlasticPoints();

	void formBb(Vector& Bb1,Vector& Bb2,Vector& Bb3);

};
#endif
