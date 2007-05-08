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

#include <Constraint.h>

int Constraint::nConstraints=0;

Constraint::Constraint()
	:DomainObject(++nConstraints),cVal(0),fTrial(0),fConvg(0)
{
}
void Constraint::setcDof(int nodeID,int dof,double coeff)
{
	static cDof newCDof;
	// Get node
	Node* pNode=pD->get<Node>(pD->getNodes(),nodeID);
	newCDof.pNode=pNode;
	// Check if the degree of freedom is activated
	///todo Throw a warning, fix this situation where no dof is constrained
	if(pNode->getActivatedDof(dof-1<0)) return;
	newCDof.dof=dof-1;
	// Add the coeffiecient
	newCDof.coeff=coeff;
	// Keep constrained dof
	theCDofs.push_back(newCDof);
}
const cDof& Constraint::getcDof(int i)
{
	return theCDofs[i];
}
void Constraint::setcVal(double val)
{
	cVal=val;
}
double Constraint::getcVal(double time)
{
	return cVal;
}
int Constraint::getncDofs()
{
	return theCDofs.size();
}
double Constraint::getDisp(int i)
{
	return theCDofs[i].pNode->getDispTrialAtDof(theCDofs[i].dof);
}
double Constraint::getVelc(int i)
{
	return theCDofs[i].pNode->getVelcTrialAtDof(theCDofs[i].dof);
}
double Constraint::getDispConvg(int i)
{
	return theCDofs[i].pNode->getDispConvgAtDof(theCDofs[i].dof);
}
double Constraint::getDispTrial(int i)
{
	return theCDofs[i].pNode->getDispTrialAtDof(theCDofs[i].dof);
}
double Constraint::getVelcConvg(int i)
{
	return theCDofs[i].pNode->getVelcConvgAtDof(theCDofs[i].dof);
}
double Constraint::getAcclConvg(int i)
{
	return theCDofs[i].pNode->getAcclConvgAtDof(theCDofs[i].dof);
}
void Constraint::update(double f)	
{
	fTrial+=f;
}
void Constraint::commit()
{
	fConvg=fTrial;
}
double Constraint::getF()
{
	return fTrial;
}
