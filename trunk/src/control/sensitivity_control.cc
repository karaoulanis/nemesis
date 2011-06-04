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

#include "control/sensitivity_control.h"
#include "analysis/analysis.h"
#include "domain/domain.h"
#include "model/model.h"
#include "model/model_element.h"
#include "model/model_node.h"
#include "soe/soe.h"

SensitivityControl::SensitivityControl()
    : ds(0),
      currParameter(0) {
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
  pA->get_soe()->zeroB();
  pA->get_domain()->zeroSensitivityParameters();
  // theLoadCase->applySensitivityParameter(currParameter);
  // Take contribution from Elements
  for (unsigned i = 0; i < pA->get_model()->get_model_elements().size(); i++) {
    ModelElement* p = pA->get_model()->get_model_elements()[i];
    this->formElementalResidual(p);
    pA->get_soe()->insertVectorIntoB(p->get_vector(), p->get_FTable(), -1.0);
  }
}
void SensitivityControl::init() {
  currParameter = 0;
  int size = pA->get_model()->get_num_eqns();
  ds.resize(size);
  ds.clear();
  for (unsigned i = 0; i < pA->get_model()->get_model_nodes().size(); i++) {
    // ModelNode* p = pA->get_model()->get_model_nodes()[i];
    /// @todo  initSensitivityMatrix stupid
    // p->get_node()->initSensitivityMatrix(theLoadCase->get_num_sens_param());
  }
}
void SensitivityControl::commit() {
  ds = pA->get_soe()->get_X();
  pA->get_model()->commitSens(ds, currParameter);
  currParameter++;
}
