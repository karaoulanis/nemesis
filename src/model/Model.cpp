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

#include <Model.h>
#include <algorithm>

/**
 * Default constructor.
 */
Model::Model(Domain* pDomain)
:AnalysisObject(),theDomain(pDomain),constrained(false)
{
}
/**
 * Destructor.
 */
Model::~Model()
{
	this->clear();
}
/**
 * Set constrained variable.
 * @param b True/false depending whether constraints have been imposed.
 */ 
void Model::setConstrained(bool b)
{
	constrained=b;
}
/**
 * Check whether constraints have been imposed.
 */ 
bool Model::isConstrained()
{
	return constrained;
}
/**
 * Set reordered variable.
 * @param b True/false depending whether reording has been done.
 */ 
void Model::setReordered(bool b)
{
	reordered=b;
}
/**
 * Check whether reording has been done.
 */ 
bool Model::isReordered()
{
	return reordered;
}
int Model::getnNodes()
{
	return theDomain->getNodes().size();
}
int Model::getnElements()
{
	return theDomain->getElements().size();
}
int Model::addModelNode(ModelNode* pModelNode)
{
	theModelNodes.push_back(pModelNode);
	return 0;
}
int Model::addModelElement(ModelElement* pModelElement)
{
	theModelElements.push_back(pModelElement);
	return 0;
}
const ModelNodeContainer& Model::getModelNodes() const
{
	return theModelNodes;
}
const ModelElementContainer& Model::getModelElements() const
{
	return theModelElements;
}
/**
 * Returns the SOE index of a dof.
 * \param NodeID The nodal id.
 * \param localDof The localDof.
 * \return The index (starting from 0) of the dof in the SOE.
 */
int Model::getSOEDof(int NodeID,int localDof)
{
	///@todo I do not like this function.
	///@bug When dof is fixed using eleimination (minus) it returnes negative.
	int n=theModelNodes.size();
	Node* pNode=0;
	ModelNode* pModelNode=0;
	for(int i=0;i<n;i++)
	{
		if(theModelNodes[i]->getNode()->getID()==NodeID)
		{
			pModelNode=getModelNodes()[i];
			pNode=pModelNode->getNode();
			break;
		}
	}
	if(pNode==0) exit(-9856);
	return pModelNode->getFTable()[pNode->getActivatedDof(localDof-1)];
}

void Model::setEquations(int n)
{
	nEquations=n;
}
int Model::getnEquations()
{		
	return nEquations;
}
void Model::incTrialDisp(const Vector& du)
{
	for(unsigned i=0;i<theModelNodes.size();i++)
		theModelNodes[i]->incTrialDisp(du);
}
void Model::setTrialDisp(const Vector& du)
{
	for(unsigned i=0;i<theModelNodes.size();i++)
		theModelNodes[i]->setTrialDisp(du);
}
void Model::incTrialVecs(const Vector& du,const Vector& da,const Vector& dv)
{
	for(unsigned i=0;i<theModelNodes.size();i++)
		theModelNodes[i]->incTrialVecs(du,da,dv);
}
void Model::setTrialVecs(const Vector& u,const Vector& a,const Vector& v)
{
	for(unsigned i=0;i<theModelNodes.size();i++)
		theModelNodes[i]->setTrialVecs(u,a,v);
}
void Model::update()
{
	for(unsigned i=0;i<theModelElements.size();i++) theModelElements[i]->update();
}
void Model::commit()
{
	for(unsigned i=0;i<theModelNodes.size();i++)	theModelNodes[i]->commit();
	for(unsigned i=0;i<theModelElements.size();i++) theModelElements[i]->commit();
}
void Model::commitSens(const Vector& ds,int param)
{
	for(unsigned i=0;i<theModelNodes.size();i++) theModelNodes[i]->commitSens(ds,param);
}
void Model::clear()
{	
	Containers::vector_delete(theModelNodes);
	Containers::vector_delete(theModelElements);
	nEquations=0;
	constrained=false;
	reordered=false;
}
int Model::getDirectedGraph(DirectedGraph& G)
{
	typedef graph_traits<DirectedGraph>::vertex_descriptor Vertex;
	typedef graph_traits<DirectedGraph>::vertices_size_type size_type;
	typedef std::pair<int,int> Pair;
	std::set<Pair> theEdges;
	Pair Edge;
	for(unsigned k=0;k<theModelElements.size();k++) 
	{
		ModelElement* pModelElem=theModelElements[k];
		for(unsigned i=0;i<pModelElem->getFTable().size();i++)
			for(unsigned j=0;j<pModelElem->getFTable().size();j++)
			{
				if(i==j) continue;
				Edge.first=pModelElem->getFTable()[i];
				Edge.second=pModelElem->getFTable()[j];
				if(Edge.first<0||Edge.second<0) continue;
				theEdges.insert(Edge);
			}
	}
	std::set<Pair>::iterator iEdges;
	for(iEdges=theEdges.begin();iEdges!=theEdges.end();iEdges++)
		add_edge(iEdges->first,iEdges->second,G);
	return 0;
}
int Model::getUndirectedGraph(UndirectedGraph& G)
{
	typedef graph_traits<UndirectedGraph>::vertex_descriptor Vertex;
	typedef graph_traits<UndirectedGraph>::vertices_size_type size_type;
	typedef std::pair<int,int> Pair;
	std::set<Pair> theEdges;
	Pair Edge;
	for(unsigned k=0;k<theModelElements.size();k++) 
	{
		ModelElement* pModelElem=theModelElements[k];
		for(unsigned i=0;i<pModelElem->getFTable().size();i++)
			for(unsigned j=i+1;j<pModelElem->getFTable().size();j++)
			{
				Edge.first=pModelElem->getFTable()[i];
				Edge.second=pModelElem->getFTable()[j];
				if(Edge.first<0||Edge.second<0) continue;
				if(Edge.first>Edge.second)std::swap(Edge.first,Edge.second);
				theEdges.insert(Edge);
			}
	}
	std::set<Pair>::iterator iEdges;
	for(iEdges=theEdges.begin();iEdges!=theEdges.end();iEdges++)
		add_edge(iEdges->first,iEdges->second,G);
	return 0;
}
void Model::setNodalStress()
{
	theDomain->zeroNodalStress();
	for(ElementIterator eIter=theDomain->getElements().begin();
		eIter!=theDomain->getElements().end();eIter++) 
			eIter->second->recoverStresses();
}
void Model::enrich()
{
	// First enrich nodes (Level set initialization)
	for(NodeIterator nIter=theDomain->getNodes().begin();
		nIter!=theDomain->getNodes().end();nIter++) 
			nIter->second->evalLevelSets();
	// And then enrich elements
	for(ElementIterator eIter=theDomain->getElements().begin();
		eIter!=theDomain->getElements().end();eIter++) 
			eIter->second->enrich();
}
