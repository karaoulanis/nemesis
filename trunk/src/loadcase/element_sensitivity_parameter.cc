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

#include "loadcase/element_sensitivity_parameter.h"

int ElementSensitivityParameter::nSensitivityParameters=0;

ElementSensitivityParameter::ElementSensitivityParameter()
	:DomainObject()
{
}
ElementSensitivityParameter::ElementSensitivityParameter(int elemID,int param)
	:DomainObject(++nSensitivityParameters)
{
	myElement=pD->get<Element>(pD->getElements(),elemID);	
	parameter=param;
}
int ElementSensitivityParameter::apply()
{
	myElement->activateParameter(parameter);
	return 0;
}