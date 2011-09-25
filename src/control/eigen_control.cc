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

#include "control/eigen_control.h"
#include "analysis/analysis.h"
#include "domain/domain.h"
#include "model/model.h"
#include "model/model_element.h"
#include "model/model_node.h"
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
void EigenControl::FormElementalTangent(ModelElement* pModelElement)  {
  pModelElement->zeroMatrix();
  pModelElement->add_K(1.0);
}

/**
 * Forms the elemental tangent for a static analysis, element by element.
 * @param pModelElement A pointer to the ModelElement that is treated.
 * @return 0 if everything is ok.
 */
void EigenControl::FormElementalMassMatrix(ModelElement* pModelElement)  {
  pModelElement->zeroMatrix();
  pModelElement->add_M(1.0);
}

/**
 * Forms the elemental residual for a static analysis, element by element.
 * @param pModelElement A pointer to the ModelElement that is treated.
 * @param time Current(?) time.
 * @return 0 if everything is ok.
 */
void EigenControl::FormElementalResidual(ModelElement* pModelElement,
                                         double /*time*/) {
  pModelElement->zeroVector();
  pModelElement->add_R(1.0);
}

/**
 * Forms the nodal tangent for a static analysis, node by node.
 * @param pModelNode A pointer to the ModelNode that is treated.
 * @return 0 if everything is ok.
 */
void EigenControl::FormNodalResidual(ModelNode* pModelNode) {
  pModelNode->zeroVector();
  pModelNode->add_R(1.0);
}

/**
 *
 */
void EigenControl::FormTangent() {
  EigenSOE* pSOE = static_cast<EigenSOE*>(pA->get_soe());
  pSOE->zeroA();
  pSOE->zeroM();

  int n = pA->get_model()->get_model_elements().size();
  for (int i = 0; i < n; i++) {
    ModelElement* p = pA->get_model()->get_model_elements()[i];
    this->FormElementalTangent(p);
    pSOE->insertMatrixIntoA(p->get_matrix(), p->get_FTable(), 1.0);
    this->FormElementalMassMatrix(p);
    pSOE->insertMatrixIntoM(p->get_matrix(), p->get_FTable(), 1.0);
  }
}

/**
 * Aborts step.
 * After a unsuccessfull step the displacement in the domain are set to the
 * previous converged step. This happens throught the ModelNodes.
 * @return 0 if everything is ok.
 */
void EigenControl::Rollback() {
}

