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

#ifndef _SPRING_H
#define _SPRING_H

#include <Element.h>
#include <SpringMaterial.h>

class Spring: public Element
{
protected:
	int nDim;
	//double gap;
	SpringMaterial* mySpringMaterial;
	Matrix T;
public:
	// Constructors and Destructor
	Spring();
	Spring(int ID,int Node_1,int Node_2,int matID,
		   double xp1=1,double xp2=0,double xp3=0,
		   double yp1=0,double yp2=1,double yp3=0);
	~Spring();
	void update();
	void commit();
	bool checkIfAllows(FEObject* f);
	void recoverStresses();
	
    const Matrix& getM();
    const Matrix& getK();
	const Vector& getR();
	const Vector& getReff();
	const Vector& getRgrad();

	// Tracker member functions
	void addTracker(int index);
	Tracker* getTracker(int index);
	void track();
};
#endif
