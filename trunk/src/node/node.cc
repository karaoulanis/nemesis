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
#include <sstream>
#include <string>
#include "domain/domain.h"
#include "tracker/tracker.h"

/**
 * Defuault Constructor.
 */
Node::Node() {
}
/**
 * Constructor.
 */
Node::Node(int ID, double xc1, double xc2, double xc3)
:DomainObject(ID),
  myActivatedDofs(MAX_NUMBER_OF_DOFS, -1),
  myConstrainedDofs(MAX_NUMBER_OF_DOFS, 0) {
  myTag = TAG_NODE;
  x1 = xc1;
  x2 = xc2;
  x3 = xc3;
  nActivatedDofs = 0;
  stress.resize(6, 0.);
  strain.resize(6, 0.);
  avgStress = 0;
  myTracker = 0;
  active_ = 0;
}
Node::~Node() {
  if (myTracker != 0) delete myTracker;
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

void Node::SetActive(bool active) {
  active ? ++active_ : --active_;
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
    velcConvg.resize(nActivatedDofs);
    velcConvg.clear();
    acclTrial.resize(nActivatedDofs);
    acclTrial.clear();
    velcTrial.resize(nActivatedDofs);
    velcTrial.clear();
    acclConvg.resize(nActivatedDofs);
    acclConvg.clear();
    dispTrial.resize(nActivatedDofs);
    dispTrial.clear();
    dispConvg.resize(nActivatedDofs);
    dispConvg.clear();
    P.resize(nActivatedDofs);
    P.clear();
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
  P.clear();
  isLoadApplied = false;
}

void Node::multDisp(double facD) {
  dispTrial*=facD;
  dispConvg*=facD;
}
void Node::zeroStress() {
  stress.clear();
  avgStress = 0;
}
void Node::addStress(const Vector& s) {
  stress+=s;
  avgStress+=1;
}
const Packet& Node::get_packet() {
  thePacket.zero();
  thePacket.tag = this->get_tag();
  thePacket.id = this->get_id();
  thePacket.dblArray[0]=x1;
  thePacket.dblArray[1]=x2;
  thePacket.dblArray[2]=x3;
  for (unsigned i = 0;i < myActivatedDofs.size();i++) {
    if (myActivatedDofs[i] >= 0) {
      thePacket.dblArray[i+3]=dispConvg[myActivatedDofs[i]];
      thePacket.dblArray[i+3+MAX_NUMBER_OF_DOFS]=velcConvg[myActivatedDofs[i]];
      thePacket.dblArray[i+3+2*MAX_NUMBER_OF_DOFS]=acclConvg[myActivatedDofs[i]];
    }
  }
  if (avgStress>0) stress*=1.0/avgStress;
  for (int i = 0;i < 6;i++)
      thePacket.dblArray[i+3+3*MAX_NUMBER_OF_DOFS]=stress[i];
  return thePacket;
}
void Node::set_packet(const Packet& p) {
  for (int i = 0; i < MAX_NUMBER_OF_DOFS; i++) {
    if (myActivatedDofs[i] < 0) continue;
    dispConvg[myActivatedDofs[i]]=p.dblArray[i+3];
    velcConvg[myActivatedDofs[i]]=p.dblArray[i+3+MAX_NUMBER_OF_DOFS];
    acclConvg[myActivatedDofs[i]]=p.dblArray[i+3+2*MAX_NUMBER_OF_DOFS];
  }
  dispTrial = dispConvg;
  velcTrial = velcConvg;
  acclTrial = acclConvg;
}
void Node::save(std::ostream& s) {
  s << "{";
  s << "\"id\":" << myID <<",";
  s << "\"crds\":[" << x1 << "," << x2 << "," << x3 << "],";
  s << "\"dofs\":[";
  /// @todo make this more general
  for (unsigned i = 0; i < myActivatedDofs.size(); i++) {
      if (i > 0) s << ',';
      s << myActivatedDofs[i];
  }
  s << "],";
  s << "\"disp\":[" << dispConvg  << "],";
  s << "\"velc\":[" << velcConvg  << "],";
  s << "\"accl\":[" << acclConvg  << "],";

  if (avgStress > 0) stress*=1.0/avgStress;
  avgStress = 0;
  s << "\"strs\":[" << stress     << "]";

  /// @todo include strains, sensitivities and eigenvalues also
  // s << "\"strn\":[" << strain     << "]";
  // s << "\"sens\":[" << dispSensi  << "],";
  // s << "\"eige\":[" << eigenVecs  << "]}";
  s << "}";
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
void Node::commit() {
  dispConvg = dispTrial;
  velcConvg = velcTrial;
  acclConvg = acclTrial;
  this->track();
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
/**
 * Add a Tracker to Node.
 * \a myTracker pointer should be \a null up to this point. If not this means
 * that a Tracker is already added and nothing is changed.
 * The Node deconstructor should take the responsibility to delete the Tracker.
 * Then track is called, in order to initialize the tracker.
 * @todo Check tracker initialization.
 */
void Node::addTracker() {
  if (myTracker != 0) return;
    myTracker = new Tracker();
  this->track();
}
/**
 * Get the Tracker.
 * An exception is thrown if no tracker is set.
 * @todo Change this to a constant pointer.
 * @return  pointer to \a myTracker.
 */
Tracker* Node::get_tracker() {
  if (myTracker == 0)
    throw SException("[nemesis:%d] No tracker is set for Node %d.", 9999, myID);
  return myTracker;
}
/**
 * Add a record to the tracker.
 * If \a myTracker pointer is null (no tracker is added) just return.
 * Otherwise gather info and send them to the tracker.
 * The domain should be already updated!
 */
void Node::track() {
  if (myTracker == 0) return;
  // define an output string stream
  std::ostringstream s;
  // start saving
  s << "{";
  // save lambda
  s << "\"lambda\":" << pD->get_lambda() << ",";
  // save time
  s << "\"time\":" << pD->get_time_curr() << ",";
  // save self
  s << "\"data\":";
  this->save(s);
  // finalize
  s << "}";
  // convert to c style string and return
  // needs to be converted to a static string before
  /// @todo: check for refactoring
  static string tmp;
  tmp = s.str();
  myTracker->track(tmp.c_str());
}
// Sensitivity functions
void Node::initSensitivityMatrix(int nGrads) {
  dispSensi.Resize(nActivatedDofs, nGrads, 0.);
}
void Node::commitSens(const Vector& v, int param) {
  dispSensi.appendCol(v, 0, param);
}
// Enrichment functions
void Node::evalLevelSets() {
}
