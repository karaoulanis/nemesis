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

#include "imposer/EliminationImposer.h"

EliminationImposer::EliminationImposer()
	:Imposer()
{
	myTag=TAG_IMPOSER_ELIMINATION;
}
EliminationImposer::~EliminationImposer()
{
}
int EliminationImposer::impose()
{
	if(theModel->isConstrained())	return 0;

	Node* pNode;									
	Element* pElement;								
	Constraint* pConstraint;
	StandardModelNode* pStdModelNode;				
	StandardModelElement* pStdModelElement;			
	EliminationModelNode* pElimModelNode;			
	EliminationModelElement* pElimModelElement;		

	//-------------------------------------------------------------------------
	// Find nodal global numbering and store it
	//-------------------------------------------------------------------------
	theModel->clear();
	int nDofs=this->createGlobalDofNumbering();

	//-------------------------------------------------------------------------
	// Number constraints with global numbering
	//-------------------------------------------------------------------------
	theNewDofs.resize(nDofs);
	for(int i=0;i<nDofs;i++) theNewDofs[i]=1;
	for(ConstraintIterator cIter=theConstraints->begin();
		cIter!=theConstraints->end();cIter++)
	{	
		pConstraint=cIter->second;;
		if(pConstraint->getncDofs()==0) continue;
		if(pConstraint->getncDofs()>1) 
			throw SException("[nemesis:%d] %s",9999,
			"Elimination method cannot be used for MultiFreedom constraints.");
		int globalDof=this->getGlobalDof(
			pConstraint->getcDof(0).pNode->getID(),pConstraint->getcDof(0).dof);
		// Check if dof is activated
		if(globalDof<0) continue;
		theNewDofs[globalDof]=-pConstraint->getID();
		if(fabs(pConstraint->getcVal())>1e-18)
			throw SException("[nemesis:%d] %s",9999,
			"Elimination method cannot be used for non homogeneous constraints.");
	}
	int nEquations=0;
	for(int i=0;i<nDofs;i++) if(theNewDofs[i]>0) theNewDofs[i]=nEquations++;
	theModel->setEquations(nEquations);

	//=========================================================================
	// Create ModelNodes
	//=========================================================================
	for(NodeIterator nIter=pA->getDomain()->getNodes().begin();
		nIter!=pA->getDomain()->getNodes().end();nIter++) 
	{
		static IDContainer oldFTable;
		static IDContainer newFTable;
		// Get a pointer to a node
		pNode=nIter->second;
		oldFTable=this->getGlobalDofs(pNode->getID());
		newFTable.resize(oldFTable.size());
		// Map old dofs to new ones
		for(unsigned i=0;i<oldFTable.size();i++) 
			newFTable[i]=theNewDofs[oldFTable[i]];
		// If constraints are found within the FTable
		if(Containers::all_positive(newFTable)) 
		{
			// Create StandardModelNode
			pStdModelNode=new StandardModelNode(newFTable,pNode);
			theModel->addModelNode(pStdModelNode);
		}
		else
		{
			// Create EliminationModelNode
			pElimModelNode=new EliminationModelNode(newFTable,pNode);
			pElimModelNode->setTheOldFTable(oldFTable);
			theModel->addModelNode(pElimModelNode);
		}
	} 
	//=========================================================================
	// Create ModelElements
	//=========================================================================
	for(ElementIterator eIter=pA->getDomain()->getElements().begin();
		eIter!=pA->getDomain()->getElements().end();eIter++) 
	{
		///@todo Needs speed improvements
		// Get next (randomly chosen) element
		pElement=eIter->second;
		// Get the ids of the nodes
		IDContainer myNodalIDs=pElement->getNodalIDs();
		// Get the local dofs of each node of the element
		IDContainer theNodalLocalDofs=pElement->getLocalNodalDofs();
		// Create an element freedom table
		IDContainer elemFTable(myNodalIDs.size()*theNodalLocalDofs.size());

		for(unsigned j=0;j<myNodalIDs.size();j++)
		{
			int NodeID=myNodalIDs[j];
			for(unsigned k=0;k<theNodalLocalDofs.size();k++)
			{
				int localDof=theNodalLocalDofs[k];
				int globalDof=this->getGlobalDof(NodeID,localDof);
				int newGlobalDof=theNewDofs[globalDof];
				elemFTable[j*theNodalLocalDofs.size()+k]=newGlobalDof;
			}
		}
		// Check if constraint found
		if(Containers::all_positive(elemFTable)) 
		{
			// Create StandardModelElement
			pStdModelElement=new StandardModelElement(elemFTable,pElement);
			theModel->addModelElement(pStdModelElement);
		}
		else
		{
			// Create EliminationModelElement
			pElimModelElement=new EliminationModelElement(elemFTable,pElement);
			theModel->addModelElement(pElimModelElement);
		}
	}
	theModel->setConstrained(true);
	return 0;
}
