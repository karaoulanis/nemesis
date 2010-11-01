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

// Included files
#include "constraints/constraint.h"

// Initialize variables
int Constraint::nConstraints=0;

/**
 * Constructor.
 */
Constraint::Constraint()
	:DomainObject(++nConstraints),cVal(0),fTrial(0),fConvg(0)
{
}
/**
 * Create a constrained dof.
 * Check if \a nodeID and \a dof is valid; otherwise throw an exception.
 * @param nodeID The nodal id to be constrained.
 * @param dof The nodal dof to be constrained.
 * @param coeff The Constraint's coefficient.
 */
void Constraint::setcDof(int nodeID,int dof,double coeff)
{
	static cDof newCDof;
	// Get node; if does not exists re-throw the exception
	Node* pNode;
	try
	{
		pNode=pD->get<Node>(pD->getNodes(),nodeID);
	}
	catch(SException)
	{
		throw;
	}
	// Check if dof is active
	if(pNode->getActivatedDof(dof-1)<0)
	{
		throw SException("[nemesis:%d] Node %d dof %d is not yet active.\n",1110,pNode->getID(),dof);
		//Containers::vector_print(pNode->getConnectedElements());cout<<endl;
		//return;
	} 
	// Now it is ok to continue
	newCDof.pNode=pNode;
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
double Constraint::getcVal(double /*time*/)
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
void Constraint::incTrialForce(double f)
{
	fTrial+=f;
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
