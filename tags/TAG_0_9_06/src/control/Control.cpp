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

#include <Control.h>
#include <cmath>

Control::Control()
{
	lambdaTrial=0;
	lambdaConvg=0;
	DLambda=0;
	dLambda=0;
}
Control::~Control()
{
}
void Control::formTangent()
{
	pA->getSOE()->zeroA();
	int n=pA->getModel()->getModelElements().size();
	for(int i=0;i<n;i++)
	{
		ModelElement* p=pA->getModel()->getModelElements()[i];
		this->formElementalTangent(p);
		pA->getSOE()->insertMatrixIntoA(p->getMatrix(),p->getFTable(),1.0);
	}
}
double Control::getLambda()
{
	return lambdaTrial;
}