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

#ifndef _TRANSIENTCONTROL_H
#define _TRANSIENTCONTROL_H

#include "control/Control.h"

class TransientControl :public Control
{
protected:
	Vector u,v,a,ut,vt,at;
	double c[3];
public:

	// Constructors and destructor
	TransientControl();
	virtual ~TransientControl();

	// Form tangent and residuals
	virtual void formElementalTangent(ModelElement* pModelElement);
	virtual void formElementalResidual(ModelElement* pModelElement,double time=0.);
	virtual void formNodalResidual(ModelNode* pModelNode);
	virtual void formResidual(double factor);
	
	virtual void init();
	virtual void commit();
	virtual void rollback();
};
#endif
