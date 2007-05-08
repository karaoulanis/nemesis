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

#include <ModelNode.h>
/**
 * Default constructor.
 */
ModelNode::ModelNode()
{
}
/**
 * Constructor.
 * Initializes the ModelObject, which in turn initializes the FEObject, passes
 * the FTable to the ModelObject and copies the address of it's node.
 */
ModelNode::ModelNode(const IDContainer& FTable,Node* pNode,Constraint* pConstraint)
	:ModelObject(FTable),myNode(pNode),myConstraint(pConstraint)
{
	// Set a pointer to the corresponding static Vector
	myVector=theStaticVectors[FTable.size()];
}
ModelNode::~ModelNode()
{
}
void ModelNode::add_uTrial(double factor)
{
	if(myNode!=0) myVector->add_cV(factor,myNode->getDispTrial());
}
void ModelNode::add_vTrial(double factor)
{
	if(myNode!=0) myVector->add_cV(factor,myNode->getVelcTrial());
}
void ModelNode::rollback()
{
	if(myNode!=0) myNode->rollback();
}
