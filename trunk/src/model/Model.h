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

#ifndef _MODEL_H
#define _MODEL_H

#include <AnalysisObject.h>
#include <Domain.h>
#include <ModelElement.h>
#include <StandardModelElement.h>
#include <ModelNode.h>
#include <ElementalLoad.h>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/profile.hpp>
#include <boost/graph/wavefront.hpp>
#include <boost/graph/graphviz.hpp>

// Using from boost nameespace
using boost::adjacency_list;
using boost::directedS;
using boost::undirectedS;
using boost::vecS;
using boost::vertex_color;
using boost::vertex_color_t;
using boost::default_color_type;
using boost::vertex_degree_t;
using boost::vertex_priority;
using boost::vertex_priority_t;
using boost::vertex_index;
using boost::vertex_index_t;
using boost::property;
using boost::property_map;
using boost::graph_traits;

typedef adjacency_list<	vecS,vecS,directedS, 
	property<vertex_color_t,default_color_type,	
	property<vertex_degree_t,int,
	property<vertex_priority_t,double > > > > DirectedGraph;
typedef adjacency_list<	vecS,vecS,undirectedS, 
	property<vertex_color_t, default_color_type,	
	property<vertex_degree_t,int,
	property<vertex_priority_t,double > > > > UndirectedGraph;

typedef std::vector<ModelNode*>						ModelNodeContainer;
typedef std::vector<ModelElement*>					ModelElementContainer;
typedef std::vector<ModelNode*>::const_iterator		ModelNodeIterator;

class Model: public AnalysisObject
{
private:
	Domain* theDomain;
	ModelNodeContainer	  theModelNodes;
	ModelElementContainer theModelElements;
	int nEquations;
	bool constrained;
	bool reordered;
public:
	// Contructor and destructor
	Model(Domain* pDomain);
	virtual ~Model();

	void setConstrained(bool b);
	bool isConstrained();
	void setReordered(bool b);
	bool isReordered();

	int getSOEDof(int NodeID,int localDof);

	// Retrieve data needed from the Domain concering Nodes and Elements
	int getnNodes();
	int getnElements();
	
	// Retrieve data from the Model concering ModelNodes and ModelElements
	int addModelNode(ModelNode* pModelNode);
	int addModelElement(ModelElement* pModelElement);
	const ModelNodeContainer& getModelNodes() const;
	const ModelElementContainer& getModelElements() const;

	// Handle the number of the system equations 
	void setEquations(int n);
	int getnEquations();

	void incTrialDisp(const Vector& du);
	void incTrialVecs(const Vector& du,const Vector& dv,const Vector& da);
	void setTrialDisp(const Vector& u);
	void setTrialVecs(const Vector& u,const Vector& v,const Vector& a);
	void update();
	void commit();
	void clear();

	void commitSens(const Vector& ds,int param);
	
	void setNodalStress();

	int getDirectedGraph(DirectedGraph& G);
	int getUndirectedGraph(UndirectedGraph& G);

	// XFem (or other type) enrichment
	void enrich();
};
#endif
