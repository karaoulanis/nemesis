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
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#include "elements/quad4.h"

Matrix Quad4::N(3, 4);
double Quad4::detJ;

Quad4::Quad4() {
}
Quad4::Quad4(int ID, int Node_1, int Node_2, int Node_3, int Node_4, int matID,
           int /*integrationRuleXi*/, int /*integrationRuleEta*/)
:Element(ID, matID) {
  myTag = TAG_ELEM_QUAD_4_DISP;
  // Get nodal data
  myNodalIDs.resize(4);
  myNodalIDs[0]=Node_1;
  myNodalIDs[1]=Node_2;
  myNodalIDs[2]=Node_3;
  myNodalIDs[3]=Node_4;
  // Set local nodal dofs
  myLocalNodalDofs.resize(2);
  myLocalNodalDofs[0]=0;
  myLocalNodalDofs[1]=1;
  // Handle common info
  this->handleCommonInfo();

  ///@todo User defined integration rule
  // p1 = integrationRuleXi;
  // p2 = integrationRuleEta;
  p1 = 2;
  p2 = 2;
  myMatPoints.resize(4);
  MultiaxialMaterial* pMat = static_cast<MultiaxialMaterial*>(myMaterial);
  myMatPoints[0]=new MatPoint(pMat, 1, 1, p1, p2);
  myMatPoints[1]=new MatPoint(pMat, 2, 1, p1, p2);
  myMatPoints[2]=new MatPoint(pMat, 2, 2, p1, p2);
  myMatPoints[3]=new MatPoint(pMat, 1, 2, p1, p2);
  for (unsigned i = 0; i < 4; i++) {
    this->findShapeFunctionsAt(myMatPoints[i]);
    double xG = N(0, 0)*x(0, 0)+N(0, 1)*x(1, 0)+N(0, 2)*x(2, 0)+N(0, 3)*x(3, 0);
    double yG = N(0, 0)*x(0, 1)+N(0, 1)*x(1, 1)+N(0, 2)*x(2, 1)+N(0, 3)*x(3, 1);
    myMatPoints[i]->setX(xG, yG);
  }
}
Quad4::~Quad4() {
  Containers::vector_delete(myMatPoints);
}
void Quad4::findShapeFunctionsAt(MatPoint* pMatPoint) {
  double xi= pMatPoint->get_r();
  double eta = pMatPoint->get_s();

  N(0, 0)=0.25*(1-xi)*(1-eta);           // N1
  N(0, 1)=0.25*(1+xi)*(1-eta);           // N2
  N(0, 2)=0.25*(1+xi)*(1+eta);           // N3
  N(0, 3)=0.25*(1-xi)*(1+eta);           // N4

  Matrix J(2, 2);
  J(0, 0)=0.25*(-x(0, 0)*(1-eta) +x(1, 0)*(1-eta) +x(2, 0)*(1+eta) -x(3, 0)*(1+eta));
  J(0, 1)=0.25*(-x(0, 0)*(1-xi)  -x(1, 0)*(1+xi)  +x(2, 0)*(1+xi)  +x(3, 0)*(1-xi));
  J(1, 0)=0.25*(-x(0, 1)*(1-eta) +x(1, 1)*(1-eta) +x(2, 1)*(1+eta) -x(3, 1)*(1+eta));
  J(1, 1)=0.25*(-x(0, 1)*(1-xi)  -x(1, 1)*(1+xi)  +x(2, 1)*(1+xi)  +x(3, 1)*(1-xi));

  detJ = J(0, 0)*J(1, 1)-J(0, 1)*J(1, 0);     // detJ
  if (detJ < 0.) cout << myNodes[0]->getID() << endl;
  double dxidx  = J(1, 1)/detJ;
  double detadx =-J(1, 0)/detJ;
  double dxidy  =-J(0, 1)/detJ;
  double detady = J(0, 0)/detJ;

  N(1, 0)=-0.25*(1-eta)*dxidx -0.25*(1-xi)*detadx;  // N1, 1
  N(1, 1)=+0.25*(1-eta)*dxidx -0.25*(1+xi)*detadx;  // N2, 1
  N(1, 2)=+0.25*(1+eta)*dxidx +0.25*(1+xi)*detadx;  // N3, 1
  N(1, 3)=-0.25*(1+eta)*dxidx +0.25*(1-xi)*detadx;  // N4, 1

  N(2, 0)=-0.25*(1-eta)*dxidy -0.25*(1-xi)*detady;  // N1, 2
  N(2, 1)=+0.25*(1-eta)*dxidy -0.25*(1+xi)*detady;  // N2, 2
  N(2, 2)=+0.25*(1+eta)*dxidy +0.25*(1+xi)*detady;  // N3, 2
  N(2, 3)=-0.25*(1+eta)*dxidy +0.25*(1-xi)*detady;  // N4, 2
}
void Quad4::commit() {
  for (unsigned int i = 0;i < myMatPoints.size();i++)
    myMatPoints[i]->getMaterial()->commit();
}
bool Quad4::checkIfAllows(FEObject* /*f*/) {
  return true;
}
void Quad4::addInitialStresses(InitialStresses* pInitialStresses) {
  if (myGroup->isActive()&&pInitialStresses->getGroupID() == myGroup->getID())
    for (unsigned i = 0;i < myMatPoints.size();i++)
      myMatPoints[i]->setInitialStresses(pInitialStresses);
}
void Quad4::recoverStresses() {
  ///@todo check
  if (p1 != 2||p2 != 2) return;
  double sq3 = 1.7320508075688772935274463415059;
  static Vector sigma(6);
  static Matrix E(4, 4);
  E(0, 0)=1+.5*sq3;   E(0, 1)=-.5;      E(0, 2)=1-.5*sq3;     E(0, 3)=-.5;
  E(1, 0)=-.5;        E(1, 1)=1+.5*sq3;   E(1, 2)=-.5;        E(1, 3)=1-.5*sq3;
  E(2, 0)=1-.5*sq3;   E(2, 1)=-.5;        E(2, 2)=1+.5*sq3;   E(2, 3)=-.5;
  E(3, 0)=-.5;        E(3, 1)=1-.5*sq3;   E(3, 2)=-.5;        E(3, 3)=1+.5*sq3;

  sigma.clear();
  for (unsigned i = 0; i < 4; i++) {
    sigma[0]=E(i, 0)*myMatPoints[0]->getMaterial()->getStress()[0]
        +E(i, 1)*myMatPoints[1]->getMaterial()->getStress()[0]
        +E(i, 2)*myMatPoints[2]->getMaterial()->getStress()[0]
        +E(i, 3)*myMatPoints[3]->getMaterial()->getStress()[0];
    sigma[1]=E(i, 0)*myMatPoints[0]->getMaterial()->getStress()[1]
        +E(i, 1)*myMatPoints[1]->getMaterial()->getStress()[1]
        +E(i, 2)*myMatPoints[2]->getMaterial()->getStress()[1]
        +E(i, 3)*myMatPoints[3]->getMaterial()->getStress()[1];
    sigma[2]=E(i, 0)*myMatPoints[0]->getMaterial()->getStress()[2]
        +E(i, 1)*myMatPoints[1]->getMaterial()->getStress()[2]
        +E(i, 2)*myMatPoints[2]->getMaterial()->getStress()[2]
        +E(i, 3)*myMatPoints[3]->getMaterial()->getStress()[2];
    sigma[3]=E(i, 0)*myMatPoints[0]->getMaterial()->getStress()[3]
        +E(i, 1)*myMatPoints[1]->getMaterial()->getStress()[3]
        +E(i, 2)*myMatPoints[2]->getMaterial()->getStress()[3]
        +E(i, 3)*myMatPoints[3]->getMaterial()->getStress()[3];
    myNodes[i]->addStress(sigma);
  }
}
/**
 * Add a Tracker to an Element's Material.
 * \a index is checked if is in \a myMatPoints range.
 * @param index The index to the Element's Material.
 */
void Quad4::addTracker(int index) {
  if (index < 0 || index>(int)myMatPoints.size()-1)
    throw SException("[nemesis:%d] %s", 9999, "Invalid index.\n");
  myMatPoints[index]->getMaterial()->addTracker();
}
/**
 * Get the Tracker with \a index.
 * \a index is checked if is in \a myMatPoints range.
 * An exception is thrown if no tracker is set.
 * @param index The index to the Element's Material.
 */
Tracker* Quad4::getTracker(int index) {
  if (index < 0 || index>(int)myMatPoints.size()-1)
    throw SException("[nemesis:%d] %s", 9999, "Invalid index.\n");
  if (myMatPoints[index]->getMaterial()->getTracker() == 0)
    throw SException("[nemesis:%d] No tracker is set for Element %d, index %d.", 9999, myID, index);
  return myMatPoints[index]->getMaterial()->getTracker();
}
/**
 * Add a record to the tracker.
 * For all non null trackers records are added.
 */
void Quad4::track() {
  for (unsigned i = 0;i < myMatPoints.size();i++)
    myMatPoints[i]->getMaterial()->track();
}
int Quad4::getnPlasticPoints() {
  int n = 0;
  for (unsigned int i = 0;i < myMatPoints.size();i++) {
    if (myMatPoints[i]->getMaterial()->isPlastic()) n++;
  }
  return n;
}
