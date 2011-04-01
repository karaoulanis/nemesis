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

#include "elements/element.h"
#include <cmath>
#include "domain/domain.h"
#include "group/group.h"
#include "loadcase/group_data.h"
#include "material/material.h"
#include "node/node.h"

Matrix** Element::theStaticMatrices = 0;
Vector** Element::theStaticVectors = 0;

/**
 * Default Constructor
 */
Element::Element() {
}

/**
 * Constructor
 * Creates an Element with a given ID and a Material ID.
 * If this is the first element to be created then a series of static Matrices
 * and Vectors is created, which are used in common for all elements. These
 * Matrices/Vectors store element properties as K, M, C, Fint etc in order to
 * reduce memory requirements and to increase speed.
 */
Element::Element(int ID, int matID)
:DomainObject(ID), myMaterial(0), activeParameter(0) {
  // Create static Matrices/vectors if they do not exist yet
  if (theStaticMatrices == 0) {
    theStaticMatrices = new Matrix*[64];
    for (int i = 1;i < 64;i++) theStaticMatrices[i]=new Matrix(i, i, 0.);
    theStaticVectors = new Vector*[64];
    for (int i = 1;i < 64;i++) theStaticVectors[i]=new Vector(i, 0.);
  }
  // Retrieve Material pointer
  myMaterial = pD->get<Material>(pD->get_materials(), matID);
}

Element::Element(int id, std::vector<Node*> nodes) {
  // Create static Matrices/vectors if they do not exist yet
  if (theStaticMatrices == 0) {
    theStaticMatrices = new Matrix*[64];
    for (int i = 1;i < 64;i++) theStaticMatrices[i]=new Matrix(i, i, 0.);
    theStaticVectors = new Vector*[64];
    for (int i = 1;i < 64;i++) theStaticVectors[i]=new Vector(i, 0.);
  }
  // Copy pointers to nodes
  myNodes = nodes;
}

Element::~Element() {
  if (theStaticMatrices != 0) {
    for (int i = 1;i < 64;i++) delete theStaticMatrices[i];
    for (int i = 1;i < 64;i++) delete theStaticVectors[i];
    delete[] theStaticMatrices;
    delete[] theStaticVectors;
    theStaticMatrices = 0;
    theStaticVectors = 0;
  }
}

const Matrix& Element::get_C() {
  static Matrix C;
  int nDofs = myNodalIDs.size()*myLocalNodalDofs.size();
  C.Resize(nDofs, nDofs, 0.);
  Vector RayleighFactors = pD->get_rayleigh_factors();
  if (RayleighFactors.size() == 0) {
    C.Clear();
  } else {
    C.add_cM(RayleighFactors[0], this->get_K());
    C.add_cM(RayleighFactors[1], this->get_M());
  }
  return C;
}

/**
 * @todo
 */
const Vector& Element::get_Reff() {
  myVector->clear();
  return *myVector;
}

/**
 * @todo
 */
const Vector& Element::get_Rgrad() {
  myVector->clear();
  return *myVector;
}

int Element::handleCommonInfo() {
  // Define number of nodes, number of dofs
  int nNodes = myNodalIDs.size();
  int nLocalDofs = myLocalNodalDofs.size();
  int nDofs = nNodes*nLocalDofs;
    myNodes.resize(nNodes);
  // Find own matrix and vector
  myMatrix = theStaticMatrices[nDofs];
  myVector = theStaticVectors[nDofs];
  // Retrieve the node pointers and check if they do exist
  for (int i = 0;i < nNodes;i++)
    myNodes[i]=pD->get <Node>(pD->get_nodes(), myNodalIDs[i]);
  // Get nodal coordinates
  x.Resize(nNodes, 3);
  for (int i = 0; i < nNodes; i++) {
    x(i, 0)=myNodes[i]->get_x1();
    x(i, 1)=myNodes[i]->get_x2();
    x(i, 2)=myNodes[i]->get_x3();
  }
  // Resize external load vector, self weigth vector
  P.resize(nDofs);
  G.resize(nDofs);
  // Inform the nodes that the corresponding Dof's must be activated
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nLocalDofs;j++)
      myNodes[i]->addDofToNode(myLocalNodalDofs[j]);
  // Create load vector
  for (int i = 0;i < nDofs;i++) P[i]=0;
  // Get Material related info
  if (myMaterial->get_id()>0) {
    // Check if the material can be assigned
    this->checkIfAllows(myMaterial);
    // Self weight
    b.resize(3);
    double g = pD->get_gravity_accl();
    b[0]=g*(pD->get_gravity_vect()[0])*(myMaterial->get_rho());
    b[1]=g*(pD->get_gravity_vect()[1])*(myMaterial->get_rho());
    b[2]=g*(pD->get_gravity_vect()[2])*(myMaterial->get_rho());
  }
  return 0;
}

/**
 * Returns the id's of the Nodes of an element.
 * @return An IDContainer of the Nodes.
 */
const IDContainer& Element::get_nodal_ids() const {
  return myNodalIDs;
}

const IDContainer& Element::get_local_nodal_dofs() const {
  return myLocalNodalDofs;
}

const std::vector<Node*>& Element::get_nodes() const {
  return myNodes;
}

const Vector& Element::get_disp_trial() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& disp=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      disp[i*nDofs+j]=myNodes[i]->get_disp_trial_at_dof(myLocalNodalDofs[j]);
  return disp;
}

const Vector& Element::get_velc_trial() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& velc=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      velc[i*nDofs+j]=myNodes[i]->get_velc_trial_at_dof(myLocalNodalDofs[j]);
  return velc;
}

const Vector& Element::get_accl_trial() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& accl=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      accl[i*nDofs+j]=myNodes[i]->get_accl_trial_at_dof(myLocalNodalDofs[j]);
  return accl;
}

const Vector& Element::get_disp_convg() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& disp=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      disp[i*nDofs+j]=myNodes[i]->get_disp_convg_at_dof(myLocalNodalDofs[j]);
  return disp;
}

const Vector& Element::get_velc_convg() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& velc=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      velc[i*nDofs+j]=myNodes[i]->get_velc_convg_at_dof(myLocalNodalDofs[j]);
  return velc;
}

const Vector& Element::get_accl_convg() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& disp=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      disp[i*nDofs+j]=myNodes[i]->get_accl_convg_at_dof(myLocalNodalDofs[j]);
  return disp;
}

/// @todo void get_disp_incrm(u)
const Vector& Element::get_disp_incrm() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& disp=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      disp[i*nDofs+j]=myNodes[i]->get_disp_trial_at_dof(myLocalNodalDofs[j])-
                      myNodes[i]->get_disp_convg_at_dof(myLocalNodalDofs[j]);
  return disp;
}

void Element::zeroLoad() {
  P.clear();
}

void Element::addLoad(const Vector& val, double fac) {
  P.add_cV(fac, val);
}

/**
 * Initial stresses.
 * They are not used by all elements. So it is simply defined here and needs to
 * be overwritten when it should be used.
 */
void Element::AddInitialStresses(int /*direction*/,
                                 double /*h1*/, double /*s1*/,
                                 double /*h2*/, double /*s2*/, double /*K0*/) {
}

void Element::addGroundMotion(int dof, double val) {
  this->get_M();
  Matrix& M=*myMatrix;
  int pos = Containers::index_find(myLocalNodalDofs, dof);
  if (pos < 0) return;
  int nNodes = myNodalIDs.size();
  int nLocalDofs = myLocalNodalDofs.size();
  int nDofs = nLocalDofs*nNodes;
    for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      P[i*nLocalDofs+pos]-=M(j, i*nLocalDofs)*val;
}

/**
 * Set GroupData
 */
void Element::SetGroupData(GroupData* groupdata) {
  groupdata_ = groupdata;
  for (unsigned i = 0; i < myNodes.size(); i++) {
    myNodes[i]->SetActive(groupdata_->active_);
  }
}

/**
 * Check if active.
 */
bool Element::IsActive() {
  return groupdata_->active_;
}

int Element::get_num_plastic_points() {
  int n = 0;
//  for (unsigned i = 0;i < theMaterialItems.size();i++)
//    if (theMaterialItems[i]->plastified()) n++;
  return n;
}

/**
 * Serialize element.
 * The following serialization follows the json format.
 * @param s An output stream (usually a stringstream)
 */
void Element::save(std::ostream& s) {
//  s << "\""                     << myID<< "\":{";
  s << "{";
  s << "\"id\":" << myID <<",";
  s << "\"tag\":"               << this->get_tag()      << ",";
  s << "\"material\":"          << myMaterial->get_id() << ",";
  s << "\"nodes\":[";
  for (unsigned i = 0;i < myNodalIDs.size(); i++) {
    if (i>0) {
      s << ',';
    }
    s << myNodalIDs[i];
  }
  s << "],";
  s << "\"inelastic_points\":"  << this->get_num_plastic_points();
  s << "}";
}

/**
 * Add a Tracker to an Element's Material.
 * This should be overwitten depending on how an Element treats its Materials
 * @param index The index to the Element's Material
 */
void Element::addTracker(int /*index*/) {
}

/**
 * Get the Tracker with \a index.
 * This should be overwitten depending on how an Element treats its Materials
 * @param index The index to the Element's Material
 */
Tracker* Element::get_tracker(int /*index*/) {
  return 0;
}

/**
 * Add a record to the tracker.
 * This should be overwitten depending on how an Element treats its Materials
 */
void Element::track() {
}

// Enrichment functions
void Element::enrich() {
}
