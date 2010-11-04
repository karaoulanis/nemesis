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

#include "control/sensitivity_control.h"

SensitivityControl::SensitivityControl() {
}
SensitivityControl::~SensitivityControl() {
}

// Form tangent and residual element by element
void SensitivityControl::formElementalTangent(ModelElement* pModelElement) {
  pModelElement->zeroMatrix();
  pModelElement->add_K(1.0);
}
void SensitivityControl::formElementalResidual(ModelElement* pModelElement,
                                               double /*time*/) {
  pModelElement->zeroVector();
  pModelElement->add_Kgrad(1.0);
}
void SensitivityControl::formResidual(double /*factor*/) {
  pA->getSOE()->zeroB();
  pA->getDomain()->zeroSensitivityParameters();
  // theLoadCase->applySensitivityParameter(currParameter);
  // Take contribution from Elements
  for (unsigned i = 0; i < pA->getModel()->getModelElements().size(); i++) {
    ModelElement* p = pA->getModel()->getModelElements()[i];
    this->formElementalResidual(p);
    pA->getSOE()->insertVectorIntoB(p->getVector(), p->getFTable(), -1.0);
  }
}
void SensitivityControl::init() {
  currParameter = 0;
  int size = pA->getModel()->getnEquations();
  ds.resize(size);
  ds.clear();
  for (unsigned i = 0; i < pA->getModel()->getModelNodes().size(); i++) {
    // ModelNode* p = pA->getModel()->getModelNodes()[i];
    ///@todo  initSensitivityMatrix stupid
    // p->getNode()->initSensitivityMatrix(theLoadCase->getnSensitivityParameters());
  }
}
void SensitivityControl::commit() {
  ds = pA->getSOE()->getX();
  pA->getModel()->commitSens(ds, currParameter);
  currParameter++;
}
