/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]     *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#ifndef NEMESIS_ELEMENTS_BRICK8_H_
#define NEMESIS_ELEMENTS_BRICK8_H_

#include "elements/element.h"
#include "material/matpoint.h"

class Brick8: public Element
{
protected:
	std::vector<MatPoint*> myMatPoints;
	static double shp[8][4][8];
	static double detJ[8];
	static std::vector<int> perm;
public:
	// Constructors and Destructor
	Brick8();
	Brick8(int ID,
				int Node_1,int Node_2,int Node_3,int Node_4,	
				int Node_5,int Node_6,int Node_7,int Node_8,
				int matID);	
	virtual ~Brick8();

	virtual const Matrix& getK();
    virtual const Matrix& getM();
	virtual const Vector& getR();

	virtual void update();
	virtual void commit();
	
	bool checkIfAllows(FEObject* f);
	void addInitialStresses(InitialStresses* pInitialStresses);	
	void recoverStresses();
	int getnPlasticPoints();

	void shapeFunctions();
	virtual void getB(Matrix& B,int node,int gPoint)=0;
};
#endif //NEMESIS_ELEMENTS_BRICK8_H_
