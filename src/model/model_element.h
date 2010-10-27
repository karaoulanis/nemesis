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

#ifndef _MODELELEMENT_H
#define _MODELELEMENT_H

#include "model/model_object.h"

class ModelElement : public ModelObject
{
private:

protected:
	Element* myElement;
	Constraint* myConstraint;
public:
	// Constructors
	ModelElement();
	ModelElement(const IDContainer& FTable,Element* pElement,Constraint* pConstraint);

	virtual ~ModelElement();

	// Access to data members
	Element* getElement()		{return myElement;}
	Constraint* getConstraint()	{return myConstraint;}

	virtual void update();
	virtual void commit();

	virtual void add_K(double factor=1.0)=0;
	virtual void add_M(double factor=1.0)=0;
	virtual void add_C(double factor=1.0)=0;
	virtual void add_R(double factor=1.0)=0;
	virtual void add_Reff(double factor=1.0)=0;
	
	virtual void add_Kgrad(double factor=1.0)	{};
};

#endif
