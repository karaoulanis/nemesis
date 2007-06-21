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

#include <Reorderer.h>
#include <fstream>
#include <SException.h>

Reorderer::Reorderer()
{
	pA->getModel()->setReordered(false);
}
Reorderer::~Reorderer()
{
}
int Reorderer::reorder()
{
	if(pA->getModel()->isReordered()) return 0;
	std::vector<int> perm;
	// Get the permutation matrix
	if(this->getPerm(perm)>0)
	{
		// Reorder Model Nodes
		for(unsigned k=0;k<pA->getModel()->getModelNodes().size();k++) 
		{
			ModelNode* pModelNode=pA->getModel()->getModelNodes()[k];
			for(unsigned i=0;i<pModelNode->getFTable().size();i++)
				if(pModelNode->getFTable()[i]>=0)
					pModelNode->setFTable(i,perm[pModelNode->getFTable()[i]]);
		}
		// Reorder Model Elements
		for(unsigned k=0;k<pA->getModel()->getModelElements().size();k++) 
		{
			ModelElement* pModelElem=pA->getModel()->getModelElements()[k];
			for(unsigned i=0;i<pModelElem->getFTable().size();i++)
				if(pModelElem->getFTable()[i]>=0)
					pModelElem->setFTable(i,perm[pModelElem->getFTable()[i]]);
		}
		cout<<"reo: Optimization returned successfully."<<endl;
	}
	else cout<<"reo: Optimization failed."<<endl;
	pA->getModel()->setReordered(true);
	return 0;
}
