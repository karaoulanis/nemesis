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
//            C.G. Panagiotopoulos (pchr@civil.auth.gr)
//*****************************************************************************

#ifndef _SENSITIVITYCONTROL_H
#define _SENSITIVITYCONTROL_H

#include "control/control.h"

class SensitivityControl :public Control
{
protected:
	Vector ds;
	int currParameter;
public:

	// Constructors and destructor
	SensitivityControl();
	virtual ~SensitivityControl();

	// Form tangent and residual element by element 
	virtual void formElementalTangent(ModelElement* pModelElement);
	virtual void formElementalResidual(ModelElement* pModelElement,double time=0.);
	
	// Form residual node by node
	void formNodalResidual(ModelNode* /*pModelNode*/)	{}

	// Methods that are used through analysis
	virtual void init();
	virtual void predict()							{}
	virtual void correct()							{}
	virtual void commit();
	virtual void rollback()							{}

	virtual void formResidual(double factor);
};

#endif
