/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]     *
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
* along with this program.  If not, see < http://www.gnu.org/licenses/>.        *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#include < cmath>
#include < stdlib.h>
#include "node/node.h"

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
}
Node::~Node() {
  if (myTracker != 0) delete myTracker;
}
//=============================================================================
// Access to data members
//=============================================================================

double Node::getx1() {
  return x1;
}
double Node::getx2() {
  return x2;
}
double Node::getx3() {
  return x3;
}
/**
 * Adds a Node ID to theConnectedElements Container.
 * This container holds all the ID's of the elements that include this Node.
 * @param  pElement A pointer to the Elemented that is connected to this Node. 
 * @return An integer indicating if everything went ok.
 */
int Node::addEleToNode(Element *pElement) {
  myConnectedElements.push_back(pElement->getID());
  return 0;
}
/**
 * Returns thegetTheConnectedElements container.
 * This container holds all the ID's of the elements that include this Node.
 * @return A pointer to the container.
 */
const IDContainer& Node::getConnectedElements() const {
  return myConnectedElements;
}
/**
 * Activates a dof in the node.
 * \param dof The dof that should be activated.
 */
int Node::addDofToNode(int dof) {
  if (myActivatedDofs.at(dof)<0)
  {
    myActivatedDofs[dof]=nActivatedDofs++;
    velcConvg.resize(nActivatedDofs); velcConvg.clear();
    acclTrial.resize(nActivatedDofs); acclTrial.clear();
    velcTrial.resize(nActivatedDofs); velcTrial.clear();
    acclConvg.resize(nActivatedDofs); acclConvg.clear();
    dispTrial.resize(nActivatedDofs); dispTrial.clear();
    dispConvg.resize(nActivatedDofs); dispConvg.clear();
    P.resize(nActivatedDofs);     P.clear();
  }
  return 0;
}
const IDContainer& Node::getActivatedDofs() const {
  return myActivatedDofs;
}
/**
 * Returns the activated dof that corresponds to the local dof.
 * Both the activated and the local dofs start from 0 and not from 1.
 * \param localDof The localdof.
 */
int Node::getActivatedDof(int localDof) const {
  ///@todo Error if localDof>MAX_NUMBER_OF_DOFS see:Constraint for e.g.
  return myActivatedDofs[localDof];
}
int Node::getnActivatedDofs() const {
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
const Vector& Node::getR() {
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
const Packet& Node::getPacket() {
  thePacket.zero();
  thePacket.tag = this->getTag();
  thePacket.id = this->getID();
  thePacket.dblArray[0]=x1;
  thePacket.dblArray[1]=x2;
  thePacket.dblArray[2]=x3;
  for (unsigned i = 0;i < myActivatedDofs.size();i++)
      if (myActivatedDofs[i]>=0)
      {
        thePacket.dblArray[i+3]=dispConvg[myActivatedDofs[i]];
        thePacket.dblArray[i+3+MAX_NUMBER_OF_DOFS]=velcConvg[myActivatedDofs[i]];
        thePacket.dblArray[i+3+2*MAX_NUMBER_OF_DOFS]=acclConvg[myActivatedDofs[i]];
      }
  if (avgStress>0) stress*=1.0/avgStress;
  for (int i = 0;i < 6;i++)
      thePacket.dblArray[i+3+3*MAX_NUMBER_OF_DOFS]=stress[i];
  return thePacket;
}
void Node::setPacket(const Packet& p) {
  for (int i = 0; i < MAX_NUMBER_OF_DOFS; i++) {
    if (myActivatedDofs[i]<0) continue;
    dispConvg[myActivatedDofs[i]]=p.dblArray[i+3];
    velcConvg[myActivatedDofs[i]]=p.dblArray[i+3+MAX_NUMBER_OF_DOFS];
    acclConvg[myActivatedDofs[i]]=p.dblArray[i+3+2*MAX_NUMBER_OF_DOFS];
  }
  dispTrial = dispConvg;
  velcTrial = velcConvg;
  acclTrial = acclConvg;
}
void Node::save(std::ostream& s) {
  static Vector crds(3);
  crds[0]=x1;
  crds[1]=x2;
  crds[2]=x3;
  s << "NODE "  <<' ';
  s << "tag " <<1000<<' '<<myTag<<' ';
  s << "id "  <<1000<<' '<<myID<<' ';
  s << "crds "  <<' '<<crds;
  s << "disp "  <<' '<<dispConvg;
  s << "velc "  <<' '<<velcConvg;
  s << "accl "  <<' '<<acclConvg;
  s << "stress "<<' '<<stress;
  s << "strain "<<' '<<strain;
  s << "dsens  "<<' '<<dispSensi;
  s << "eigen  "<<' '<<eigenVecs;
  s << "END "<<' ';
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
void Node::setTrialDisp(const Vector& u) {
  dispTrial = u;
}
void Node::setTrialVelc(const Vector& v) {
  velcTrial = v;
}
void Node::setTrialAccl(const Vector& a) {
  acclTrial = a;
}
void Node::commit() {
  dispConvg = dispTrial;
  velcConvg = velcTrial;
  acclConvg = acclTrial;
  this->track();
}
const Vector& Node::getDispTrial() {
  return dispTrial;
}
const Vector& Node::getDispConvg() {
  return dispConvg;
}
const Vector& Node::getVelcTrial() {
  return velcTrial;
}
const Vector& Node::getVelcConvg() {
  return velcConvg;
}
const Vector& Node::getAcclTrial() {
  return acclTrial;
}
const Vector& Node::getAcclConvg() {
  return acclConvg;
}
void Node::rollback() {
  dispTrial = dispConvg;
}
double Node::getDispTrialAtDof(int dof) {
  return dispTrial[myActivatedDofs[dof]];
}
double Node::getVelcTrialAtDof(int dof) {
  return velcTrial[myActivatedDofs[dof]];
}
double Node::getAcclTrialAtDof(int dof) {
  return acclTrial[myActivatedDofs[dof]];
}
double Node::getDispConvgAtDof(int dof) {
  return dispConvg[myActivatedDofs[dof]];
}
double Node::getVelcConvgAtDof(int dof) {
  return velcConvg[myActivatedDofs[dof]];
}
double Node::getAcclConvgAtDof(int dof) {
  return acclConvg[myActivatedDofs[dof]];
}
bool Node::isActive() {
  for (unsigned i = 0;i < myConnectedElements.size();i++)
    if (pD->get < Element>(pD->getElements(),
      myConnectedElements[i])->isActive()) return true;
  return false;
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
Tracker* Node::getTracker() {
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
  ostringstream s;
  s << "DATA "  <<' ';
  s << "disp "  <<' '<<dispConvg;
  s << "velc "  <<' '<<velcConvg;
  s << "accl "  <<' '<<acclConvg;
  s << "END "<<' ';
  myTracker->track(pD->getLambda(), pD->getTimeCurr(), s.str());
}
// Sensitivity functions
void Node::initSensitivityMatrix(int nGrads) {
  dispSensi.resize(nActivatedDofs, nGrads, 0.);
}
void Node::commitSens(const Vector& v, int param) {
  dispSensi.appendCol(v, 0, param);
}
// Enrichment functions
void Node::evalLevelSets() {
}
