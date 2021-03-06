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

#include "node/node.h"
#include <stdlib.h>
#include <cmath>
#include <string>

Node::Node()
    : x1(0.),
      x2(0.),
      x3(0.),
      myConnectedElements(0),
      myActivatedDofs(MAX_NUMBER_OF_DOFS, -1),
      myConstrainedDofs(MAX_NUMBER_OF_DOFS, 0),
      nActivatedDofs(0),
      P(),
      dispTrial(),
      dispConvg(),
      velcTrial(),
      velcConvg(),
      acclTrial(),
      acclConvg(),
      dispSensi(),
      eigenVecs(),
      avgStress(0),
      stress(6, 0.),
      strain(6, 0.),
      isLoadApplied(false),
      active_(0) {
}


Node::Node(int ID, double xc1, double xc2, double xc3)
    : DomainObject(ID),
      x1(xc1),
      x2(xc2),
      x3(xc3),
      myConnectedElements(0),
      myActivatedDofs(MAX_NUMBER_OF_DOFS, -1),
      myConstrainedDofs(MAX_NUMBER_OF_DOFS, 0),
      nActivatedDofs(0),
      P(),
      dispTrial(),
      dispConvg(),
      velcTrial(),
      velcConvg(),
      acclTrial(),
      acclConvg(),
      dispSensi(),
      eigenVecs(),
      avgStress(0),
      stress(6, 0.),
      strain(6, 0.),
      isLoadApplied(false),
      active_(0) {
}


Node::~Node() {
}
//=============================================================================
// Access to data members
//=============================================================================

double Node::get_x1() {
  return x1;
}
double Node::get_x2() {
  return x2;
}
double Node::get_x3() {
  return x3;
}

void Node::set_active(bool active) {
  active_ = active;
}

void Node::AddActive(bool active) {
  active ? ++active_ : active_;
}

bool Node::IsActive() {
  return active_ >0 ? true : false;
}

/**
 * Activates a dof in the node.
 * \param dof The dof that should be activated.
 */
int Node::addDofToNode(int dof) {
  if (myActivatedDofs.at(dof) < 0) {
    /// @todo use resize(size, 0.)
    myActivatedDofs[dof]=nActivatedDofs++;
    velcConvg.Resize(nActivatedDofs);
    velcConvg.Clear();
    acclTrial.Resize(nActivatedDofs);
    acclTrial.Clear();
    velcTrial.Resize(nActivatedDofs);
    velcTrial.Clear();
    acclConvg.Resize(nActivatedDofs);
    acclConvg.Clear();
    dispTrial.Resize(nActivatedDofs);
    dispTrial.Clear();
    dispConvg.Resize(nActivatedDofs);
    dispConvg.Clear();
    P.Resize(nActivatedDofs);
    P.Clear();
  }
  return 0;
}
const IDContainer& Node::get_activated_dofs() const {
  return myActivatedDofs;
}
/**
 * Returns the activated dof that corresponds to the local dof.
 * Both the activated and the local dofs start from 0 and not from 1.
 * \param localDof The localdof.
 */
int Node::get_activated_dof(int localDof) const {
  /// @todo Error if localDof>MAX_NUMBER_OF_DOFS see:Constraint for e.g.
  return myActivatedDofs[localDof];
}
int Node::get_num_activated_dofs() const {
  return nActivatedDofs;
}
void Node::addLoad(int dof, double value, double factor) {
  P[myActivatedDofs[dof]]-=factor*value;
  isLoadApplied = true;
}
void Node::addInitialDisp(int dof, double disp) {
  dispTrial[myActivatedDofs[dof]]=disp;
}
void Node::addInitialVelc(int dof, double velc) {
  velcTrial[myActivatedDofs[dof]]=velc;
}
const Vector& Node::get_R() {
  return P;
}
bool Node::existsLoad() {
  return isLoadApplied;
}
void Node::zeroLoad() {
  P.Clear();
  isLoadApplied = false;
}

void Node::multDisp(double facD) {
  dispTrial*=facD;
  dispConvg*=facD;
}
void Node::zeroStress() {
  stress.Clear();
  avgStress = 0;
}
void Node::addStress(const Vector& s) {
  stress+=s;
  avgStress+=1;
}
void Node::AverageStresses() {
  if (avgStress > 0) {
    stress *= 1.0 / avgStress;
  }
  avgStress = 0;
}

const Vector& Node::get_stress() {
  return stress;
}

void Node::Save(std::ostream* os) {
  (*os) << "{";
  (*os) << "\"id\":" << id_ <<",";
  (*os) << "\"crds\":[" << x1 << "," << x2 << "," << x3 << "],";
  (*os) << "\"dofs\":[";
  /// @todo make this more general
  for (unsigned i = 0; i < myActivatedDofs.size(); i++) {
      if (i > 0) (*os) << ',';
      (*os) << myActivatedDofs[i];
  }
  (*os) << "],";
  (*os) << "\"disp\":[" << dispConvg  << "],";
  (*os) << "\"velc\":[" << velcConvg  << "],";
  (*os) << "\"accl\":[" << acclConvg  << "],";
  (*os) << "\"strs\":[" << stress     << "]";

  /// @todo include strains, sensitivities and eigenvalues also
  // s << "\"strn\":[" << strain     << "]";
  // s << "\"sens\":[" << dispSensi  << "],";
  // s << "\"eige\":[" << eigenVecs  << "]}";
  (*os) << "}";
}


void Node::incTrialDisp(const Vector& du) {
  dispTrial+=du;
}


void Node::addTrialVelc(const Vector& dv) {
  velcTrial+=dv;
}


void Node::addTrialAccl(const Vector& da) {
  acclTrial+=da;
}


void Node::set_trial_disp(const Vector& u) {
  dispTrial = u;
}


void Node::set_trial_velc(const Vector& v) {
  velcTrial = v;
}
void Node::set_trial_accl(const Vector& a) {
  acclTrial = a;
}
void Node::Commit() {
  dispConvg = dispTrial;
  velcConvg = velcTrial;
  acclConvg = acclTrial;
}
const Vector& Node::get_disp_trial() {
  return dispTrial;
}
const Vector& Node::get_disp_convg() {
  return dispConvg;
}
const Vector& Node::get_velc_trial() {
  return velcTrial;
}
const Vector& Node::get_velc_convg() {
  return velcConvg;
}
const Vector& Node::get_accl_trial() {
  return acclTrial;
}
const Vector& Node::get_accl_convg() {
  return acclConvg;
}
void Node::rollback() {
  dispTrial = dispConvg;
}
double Node::get_disp_trial_at_dof(int dof) {
  return dispTrial[myActivatedDofs[dof]];
}

double Node::get_velc_trial_at_dof(int dof) {
  return velcTrial[myActivatedDofs[dof]];
}

double Node::get_accl_trial_at_dof(int dof) {
  return acclTrial[myActivatedDofs[dof]];
}

double Node::get_disp_convg_at_dof(int dof) {
  return dispConvg[myActivatedDofs[dof]];
}

double Node::get_velc_convg_at_dof(int dof) {
  return velcConvg[myActivatedDofs[dof]];
}

double Node::get_accl_convg_at_dof(int dof) {
  return acclConvg[myActivatedDofs[dof]];
}

// Sensitivity functions
void Node::initSensitivityMatrix(int nGrads) {
  dispSensi.Resize(nActivatedDofs, nGrads, 0.);
}

void Node::commitSens(const Vector& v, int param) {
  dispSensi.AppendCol(v, 0, param);
}

// Enrichment functions
void Node::evalLevelSets() {
}
