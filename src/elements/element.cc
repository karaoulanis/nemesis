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
#include "group/group.h"
#include "group/group_data.h"
#include "material/material.h"
#include "node/node.h"

Matrix** Element::theStaticMatrices = 0;
Vector** Element::theStaticVectors = 0;
double Element::gravitydirection_[3] = {0., -1., 0.};
double Element::gravityacceleration_ = 9.81;
Vector Element::rayleigh_(0);


Element::Element()
    : nodes_(0) {
}

Element::Element(int id, std::vector<Node*> nodes)
    : DomainObject(id),
      nodes_(nodes) {
  // Create static Matrices/vectors if they do not exist yet
  if (theStaticMatrices == 0) {
    theStaticMatrices = new Matrix*[64];
    for (int i = 1;i < 64;i++) {
      theStaticMatrices[i]=new Matrix(i, i, 0.);
    }
    theStaticVectors = new Vector*[64];
    for (int i = 1;i < 64;i++) {
      theStaticVectors[i]=new Vector(i, 0.);
    }
  }
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
  if (rayleigh_.get_size() == 0) {
    C.Clear();
  } else {
    C.Add_cM(rayleigh_[0], this->get_K());
    C.Add_cM(rayleigh_[1], this->get_M());
  }
  return C;
}

/**
 * @todo
 */
const Vector& Element::get_Reff() {
  myVector->Clear();
  return *myVector;
}

/**
 * @todo
 */
const Vector& Element::get_Rgrad() {
  myVector->Clear();
  return *myVector;
}

void Element::set_gravitydirection(double xG, double yG, double zG) {
    gravitydirection_[0] = xG;
    gravitydirection_[1] = yG;
    gravitydirection_[2] = zG;
}

void Element::set_gravityacceleration(double acceleration) {
  gravityacceleration_ = acceleration;
}

void Element::AssignGravityLoads() {
  // Check if a Material exists
  /// @todo Provide a better check
  if (myMaterial->get_id() > 0) {
    b.Resize(3);
    double gamma = (myMaterial->get_rho())*gravityacceleration_;
    b[0] = gamma*gravitydirection_[0];
    b[1] = gamma*gravitydirection_[1];
    b[2] = gamma*gravitydirection_[2];
  }
}

void Element::set_rayleigh(const Vector& rayleigh) {
  rayleigh_.Resize(rayleigh.get_size());
  rayleigh_ = rayleigh;
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
  return nodes_;
}

const Vector& Element::get_disp_trial() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& disp=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      disp[i*nDofs+j] = nodes_[i]->get_disp_trial_at_dof(myLocalNodalDofs[j]);
  return disp;
}

const Vector& Element::get_velc_trial() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& velc=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      velc[i*nDofs+j] = nodes_[i]->get_velc_trial_at_dof(myLocalNodalDofs[j]);
  return velc;
}

const Vector& Element::get_accl_trial() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& accl=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      accl[i*nDofs+j] = nodes_[i]->get_accl_trial_at_dof(myLocalNodalDofs[j]);
  return accl;
}

const Vector& Element::get_disp_convg() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& disp=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      disp[i*nDofs+j] = nodes_[i]->get_disp_convg_at_dof(myLocalNodalDofs[j]);
  return disp;
}

const Vector& Element::get_velc_convg() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& velc=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      velc[i*nDofs+j] = nodes_[i]->get_velc_convg_at_dof(myLocalNodalDofs[j]);
  return velc;
}

const Vector& Element::get_accl_convg() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& disp=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      disp[i*nDofs+j] = nodes_[i]->get_accl_convg_at_dof(myLocalNodalDofs[j]);
  return disp;
}

/// @todo void get_disp_incrm(u)
const Vector& Element::get_disp_incrm() {
  int nDofs = myLocalNodalDofs.size();
  int nNodes = myNodalIDs.size();
  Vector& disp=*myVector;
  for (int i = 0;i < nNodes;i++)
    for (int j = 0;j < nDofs;j++)
      disp[i*nDofs+j] = nodes_[i]->get_disp_trial_at_dof(myLocalNodalDofs[j])-
                        nodes_[i]->get_disp_convg_at_dof(myLocalNodalDofs[j]);
  return disp;
}

void Element::zeroLoad() {
  P.Clear();
}

void Element::addLoad(const Vector& val, double fac) {
  P.Add_cV(fac, val);
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
void Element::SetGroupData(const GroupData* groupdata) {
  groupdata_ = groupdata;
}

/**
 * Check if active.
 */
bool Element::IsActive() {
  return groupdata_->active;
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
void Element::Save(std::ostream* os) {
  (*os) << "{";
  (*os) << "\"id\":" << id_ <<",";
  (*os) << "\"material\":"          << myMaterial->get_id() << ",";
  (*os) << "\"nodes\":[";
  for (unsigned i = 0;i < nodes_.size(); i++) {
    if (i>0) {
      (*os) << ',';
    }
    (*os) << nodes_[i]->get_id();
  }
  (*os) << "],";
  (*os) << "\"inelastic_points\":"  << this->get_num_plastic_points();
  (*os) << "}";
}

// Enrichment functions
void Element::enrich() {
}
