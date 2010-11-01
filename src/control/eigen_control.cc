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

#include "control/eigen_control.h"
#include "soe/eigen_soe.h"
/**
 * Constructor.
 */
EigenControl::EigenControl() {
}
/**
 * Destructor.
 */
EigenControl::~EigenControl() {
}
/**
 * Forms the elemental tangent for a static analysis, element by element.
 * @param pModelElement A pointer to the ModelElement that is treated.
 * @return 0 if everything is ok.
 */
void EigenControl::formElementalTangent(ModelElement* pModelElement)  {
  pModelElement->zeroMatrix();
  pModelElement->add_K(1.0);
}
/**
 * Forms the elemental tangent for a static analysis, element by element.
 * @param pModelElement A pointer to the ModelElement that is treated.
 * @return 0 if everything is ok.
 */
void EigenControl::formElementalMassMatrix(ModelElement* pModelElement)  {
  pModelElement->zeroMatrix();
  pModelElement->add_M(1.0);
}
/**
 * Forms the elemental residual for a static analysis, element by element.
 * @param pModelElement A pointer to the ModelElement that is treated.
 * @param time Current(?) time.
 * @return 0 if everything is ok.
 */
void EigenControl::formElementalResidual(ModelElement* pModelElement, double /*time*/) {
  pModelElement->zeroVector();
  pModelElement->add_R(1.0);
}
/**
 * Forms the nodal tangent for a static analysis, node by node.
 * @param pModelNode A pointer to the ModelNode that is treated.
 * @return 0 if everything is ok.
 */
void EigenControl::formNodalResidual(ModelNode* pModelNode) {
  pModelNode->zeroVector();
  pModelNode->add_R(1.0);
}
void EigenControl::formTangent() {
  EigenSOE* pSOE = static_cast < EigenSOE*>(pA->getSOE());
  pSOE->zeroA();
  pSOE->zeroM();

  int n = pA->getModel()->getModelElements().size();
  for (int i = 0; i < n; i++) {
    ModelElement* p = pA->getModel()->getModelElements()[i];
    this->formElementalTangent(p);
    pSOE->insertMatrixIntoA(p->getMatrix(), p->getFTable(), 1.0);
    this->formElementalMassMatrix(p);
    pSOE->insertMatrixIntoM(p->getMatrix(), p->getFTable(), 1.0);
  }
}
/**
 * Aborts step.
 * After a unsuccessfull step the displacement in the domain are set to the 
 * previous converged step. This happens throught the ModelNodes. 
 * @return 0 if everything is ok.
 */
void EigenControl::rollback() {
}

