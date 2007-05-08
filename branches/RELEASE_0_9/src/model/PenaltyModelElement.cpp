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


#include <PenaltyModelElement.h>

/**
 * Default constructor
 */
PenaltyModelElement::PenaltyModelElement()
	:ModelElement()
{
}
/**
 * Constructor
 */
PenaltyModelElement::PenaltyModelElement(const IDContainer& FTable,
											 Constraint* pConstraint,
											 double aFactor)
	:ModelElement(FTable,0,pConstraint)
{
	a=aFactor;
	myMatrix=theStaticMatrices[FTable.size()];
	myVector=theStaticVectors[FTable.size()];
}
/**
 * Destructor
 */
PenaltyModelElement::~PenaltyModelElement()
{
}
void PenaltyModelElement::add_K(double factor)
{
	for(unsigned int i=0;i<theFTable.size();i++)
	{
		for(unsigned int j=0;j<theFTable.size();j++)
		{
			double ci=myConstraint->getcDof(i).coeff;
			double cj=myConstraint->getcDof(j).coeff;
			(*myMatrix)(i,j)+=factor*a*ci*cj;
		}
	}
}
void PenaltyModelElement::add_M(double factor)
{
	for(unsigned int i=0;i<theFTable.size();i++)
	{
		for(unsigned int j=0;j<theFTable.size();j++)
		{
			double ci=myConstraint->getcDof(i).coeff;
			double cj=myConstraint->getcDof(j).coeff;
			(*myMatrix)(i,j)+=factor*a*ci*cj;
		}
	}
}
void PenaltyModelElement::add_C(double factor)
{
}
void PenaltyModelElement::add_R(double factor)
{
	///@todo Make this work for the non-linear multi constraint
	double c=myConstraint->getcVal(0.);
	for(unsigned i=0;i<theFTable.size();i++)
	{
		double ci=myConstraint->getcDof(i).coeff;
		double ui=myConstraint->getDisp(i);
		(*myVector)[i]=+factor*ci*(a*ui-a*c);
	}
}
void PenaltyModelElement::add_Reff(double factor)
{
	///@todo Make this work for the non-linear multi constraint
	double c=myConstraint->getcVal(0.);
	for(unsigned i=0;i<theFTable.size();i++)
	{
		double ci=myConstraint->getcDof(i).coeff;
		double ui=myConstraint->getDisp(i);
		(*myVector)[i]=+factor*ci*(a*ui-a*c);
	}
}
