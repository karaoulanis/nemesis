/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 2, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License along   *
*   with this program; if not, write to the Free Software Foundation, Inc.,   *
*   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.               *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//*****************************************************************************

#ifndef _BAR_H
#define _BAR_H

#include <Element.h>
#include <UniaxialMaterial.h>

class Bar: public Element
{
protected:
	int nDim;
	CrossSection* iSection;
	CrossSection* jSection;
	double L0;
	double A0;
	Vector cosX;
	UniaxialMaterial* myUniMaterial;
public:
	// Constructors and Destructor
	Bar();
	Bar(int ID,int Node_1,int Node_2,int matID,int iSecID,int jSecID);
	~Bar();
	void update();
	void commit();
	bool checkIfAllows(FEObject* f);
    const Matrix& getM();
	const Vector& getReff();
	void recoverStresses();
	
	// Tracker member functions
	void addTracker(int index);
	Tracker* getTracker(int index);
	void track();
};
#endif
