/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2011 F.E.Karaoulanis [http://www.nemesis-project.org]     *
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
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#include "model/model.h"
#include <algorithm>
#include <set>
#include <utility>
#include "domain/domain.h"
#include "model/model_element.h"
#include "model/model_node.h"
#include "model/standard_model_element.h"
#include "node/node.h"

/**
 * Default constructor.
 */
Model::Model(Domain* pDomain)
    : theDomain(pDomain),
      theModelNodes(0),
      theModelElements(0),
      nEquations(0),
      constrained(false),
      reordered(false) {
}
/**
 * Destructor.
 */
Model::~Model() {
  this->clear();
}
/**
 * Set constrained variable.
 * @param b True/false depending whether constraints have been imposed.
 */
void Model::set_constrained(bool b) {
  constrained = b;
}
/**
 * Check whether constraints have been imposed.
 */
bool Model::isConstrained() {
  return constrained;
}
/**
 * Set reordered variable.
 * @param b True/false depending whether reording has been done.
 */
void Model::set_reordered(bool b) {
  reordered = b;
}
/**
 * Check whether reording has been done.
 */
bool Model::isReordered() {
  return reordered;
}
int Model::get_num_nodes() {
  return theDomain->get_nodes().size();
}
int Model::get_num_elements() {
  return theDomain->get_elements().size();
}
int Model::addModelNode(ModelNode* pModelNode) {
  theModelNodes.push_back(pModelNode);
  return 0;
}
int Model::addModelElement(ModelElement* pModelElement) {
  theModelElements.push_back(pModelElement);
  return 0;
}
const ModelNodeContainer& Model::get_model_nodes() const {
  return theModelNodes;
}
const ModelElementContainer& Model::get_model_elements() const {
  return theModelElements;
}
/**
 * Returns the SOE index of a dof.
 * \param NodeID The nodal id.
 * \param localDof The localDof.
 * \return The index (starting from 0) of the dof in the SOE.
 */
int Model::get_soe_dof(int NodeID, int localDof) {
  /// @todo I do not like this function.
  /// @bug When dof is fixed using eleimination (minus) it returnes negative.
  int n = theModelNodes.size();
  Node* pNode = 0;
  ModelNode* pModelNode = 0;
  for (int i = 0; i < n; i++) {
    if (theModelNodes[i]->get_node()->get_id() == NodeID) {
      pModelNode = get_model_nodes()[i];
      pNode = pModelNode->get_node();
      break;
    }
  }
  if (pNode == 0) exit(-9856);
  return pModelNode->get_FTable()[pNode->get_activated_dof(localDof-1)];
}

void Model::set_equations(int n) {
  nEquations = n;
}
int Model::get_num_eqns() {
  return nEquations;
}
void Model::incTrialDisp(const Vector& du) {
  for (unsigned i = 0;i < theModelNodes.size();i++)
    theModelNodes[i]->incTrialDisp(du);
}
void Model::set_trial_disp(const Vector& du) {
  for (unsigned i = 0;i < theModelNodes.size();i++)
    theModelNodes[i]->set_trial_disp(du);
}
void Model::incTrialVecs(const Vector& du, const Vector& da, const Vector& dv) {
  for (unsigned i = 0;i < theModelNodes.size();i++)
    theModelNodes[i]->incTrialVecs(du, da, dv);
}
void Model::set_trial_vecs(const Vector& u, const Vector& a, const Vector& v) {
  for (unsigned i = 0;i < theModelNodes.size();i++)
    theModelNodes[i]->set_trial_vecs(u, a, v);
}
void Model::update() {
  for (unsigned i = 0;i < theModelElements.size();i++)
    theModelElements[i]->update();
}
void Model::commit() {
  for (unsigned i = 0;i < theModelNodes.size();i++)
    theModelNodes[i]->commit();
  for (unsigned i = 0;i < theModelElements.size();i++)
    theModelElements[i]->commit();
}
void Model::commitSens(const Vector& ds, int param) {
  for (unsigned i = 0;i < theModelNodes.size();i++)
    theModelNodes[i]->commitSens(ds, param);
}
void Model::clear() {
  Containers::vector_delete(theModelNodes);
  Containers::vector_delete(theModelElements);
  nEquations = 0;
  constrained = false;
  reordered = false;
}
int Model::get_directed_graph(DirectedGraph& G) {
  typedef graph_traits < DirectedGraph>::vertex_descriptor Vertex;
  typedef graph_traits < DirectedGraph>::vertices_size_type size_type;
  typedef std::pair < int, int > Pair;
  std::set < Pair > theEdges;
  Pair Edge;
  for (unsigned k = 0; k < theModelElements.size(); k++) {
    ModelElement* pModelElem = theModelElements[k];
    for (unsigned i = 0;i < pModelElem->get_FTable().size();i++)
      for (unsigned j = 0; j < pModelElem->get_FTable().size(); j++) {
        if (i == j) continue;
        Edge.first = pModelElem->get_FTable()[i];
        Edge.second = pModelElem->get_FTable()[j];
        if (Edge.first < 0||Edge.second < 0) continue;
        theEdges.insert(Edge);
      }
  }
  std::set < Pair>::iterator iEdges;
  for (iEdges = theEdges.begin();iEdges != theEdges.end();iEdges++)
    add_edge(iEdges->first, iEdges->second, G);
  return 0;
}
int Model::get_undirected_graph(UndirectedGraph& G) {
  typedef graph_traits<UndirectedGraph>::vertex_descriptor Vertex;
  typedef graph_traits<UndirectedGraph>::vertices_size_type size_type;
  typedef std::pair<int, int> Pair;
  std::set<Pair> theEdges;
  Pair Edge;
  for (unsigned k = 0; k < theModelElements.size(); k++) {
    ModelElement* pModelElem = theModelElements[k];
    for (unsigned i = 0;i < pModelElem->get_FTable().size();i++)
      for (unsigned j = i+1; j < pModelElem->get_FTable().size(); j++) {
        Edge.first = pModelElem->get_FTable()[i];
        Edge.second = pModelElem->get_FTable()[j];
        if (Edge.first < 0||Edge.second < 0) continue;
        if (Edge.first>Edge.second)std::swap(Edge.first, Edge.second);
        theEdges.insert(Edge);
      }
  }
  std::set<Pair>::iterator iEdges;
  for (iEdges = theEdges.begin(); iEdges != theEdges.end(); iEdges++)
    add_edge(iEdges->first, iEdges->second, G);
  return 0;
}
void Model::set_nodal_stress() {
  theDomain->zeroNodalStress();
  for (ElementIterator eIter = theDomain->get_elements().begin();
    eIter != theDomain->get_elements().end(); eIter++)
      eIter->second->recoverStresses();
}
void Model::enrich() {
  // First enrich nodes (Level set initialization)
  for (NodeIterator nIter = theDomain->get_nodes().begin();
    nIter != theDomain->get_nodes().end(); nIter++)
      nIter->second->evalLevelSets();
  // And then enrich elements
  for (ElementIterator eIter = theDomain->get_elements().begin();
    eIter != theDomain->get_elements().end(); eIter++)
      eIter->second->enrich();
}
