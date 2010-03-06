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

#include <Imposer.h>

Imposer::Imposer()
{
	theModel=pA->getModel();
	theConstraints=&(pA->getDomain()->getConstraints());
	theModel->setConstrained(false);
	theDomain=pA->getDomain();
}
Imposer::~Imposer()
{
}
/**
 * Creates global dof numbering and stores it into the appropriate containers.
 * There exist two containers: myNodalIDs and theNodalGlobalDofs. The first 
 * one holds the ids, while the second one holds the number of the 
 * corresponding activated dofs in an accumulative way which makes the 
 * retrieval of global dof numbering an easy task. The storing is as follows:
 * Suppose three nodes with ids and activated dofs as, (1,3), (3,2), (10,4).
 * Then it while be stored:
 * myNodalIDs        : [1 3 10]
 * theNodalGlobalDofs : [3 5  9]
 * @return The number of global dofs.
 */
int Imposer::createGlobalDofNumbering()
{
	int nNodes=theDomain->getNodes().size();
	myNodalIDs.resize(nNodes);
	theNodalGlobalDofs.resize(nNodes);
	int k=0;
	for(NodeIterator nIter=theDomain->getNodes().begin();
		nIter!=theDomain->getNodes().end();nIter++)
	{
		// Get the next node from the node container stored in the domain 
		Node* pNode=nIter->second;
		// Store the nodal id
		myNodalIDs[k]=pNode->getID();
		// Sore the number of activated dofs
		int nActivatedDofs=pNode->getnActivatedDofs();
		if(k==0)	
			theNodalGlobalDofs[0]=nActivatedDofs;
		else		
			theNodalGlobalDofs[k]=theNodalGlobalDofs[k-1]+nActivatedDofs;
		k++;
	}
	return theNodalGlobalDofs[theNodalGlobalDofs.size()-1];
}
/**
 * Retrieves a global dof when knowing the node id and the local dof id.
 * @see Imposer::createGlobalDofNumbering()
 * For this example the retrieval can be done as follows:
 * Suppose we want the global dof for node 3, local dof 2, where local dof 2 is
 * stored at the first position of the activated dofs.
 * \li Find the index of 3 in myNodalIDs: Here is 1. \n
 * \li Find starting globalDof. Here we start at theNodalGlobalDofs[1-1]=3. \n
 * \li Find position of asked dof in theActivatedDofs and add it to globalDof.
 *     Here position is 0.\n
 * \param NodeID The nodal id.
 * \param localDof The local dof that we are interested in.
 * \return -1 if dof not activated; else globalDof
 */
int Imposer::getGlobalDof(int NodeID,int localDof)
{
	int index=Containers::index_find(myNodalIDs,NodeID);
	int globalDof=0;
	if(index>0)	globalDof=theNodalGlobalDofs[index-1];
	Node* pNode=theDomain->get<Node>(theDomain->getNodes(),NodeID);
	if(pNode->getActivatedDof(localDof)<0) return -1;
	globalDof+=pNode->getActivatedDof(localDof);
	return globalDof;
}
/**
 * Retrieves a global dof when knowing the node id and the local dof id.
 * @see Model::createGlobalDofNumbering()
 * @see Model::getGlobalDof(int NodeID,int localDof)\n
 * The process is the same as mentioned above, but now all the global dofs of
 * a node are returned.
 * @param NodeID The nodal id.
 * @return A constant reference to an IDContainer holding the global dofs.
 */
const IDContainer Imposer::getGlobalDofs(int NodeID)
{
	int index=Containers::index_find(myNodalIDs,NodeID);
	int nActivatedDofs=
		theDomain->get<Node>(theDomain->getNodes(),NodeID)->getnActivatedDofs();
	IDContainer globalDofs(nActivatedDofs);
	int startDof;
	if(index==0)	startDof=0;
	else			startDof=theNodalGlobalDofs[index-1];
	for(int i=0;i<nActivatedDofs;i++) globalDofs[i]=startDof+i;
	return globalDofs;
}
