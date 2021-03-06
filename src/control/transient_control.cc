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

#include "control/transient_control.h"
#include "analysis/analysis.h"
#include "domain/domain.h"
#include "model/model_element.h"
#include "model/model_node.h"
#include "soe/soe.h"

/**
 * Constructor.
 */
TransientControl::TransientControl()
    : u(0),
      v(0),
      a(0),
      ut(0),
      vt(0),
      at(0) {
}

/**
 * Destructor.
 */
TransientControl::~TransientControl() {
}

/**
 * @return 0 if everything is ok.
 */
void TransientControl::Init() {
  int size = pA->get_model()->get_num_eqns();
  /// @todo: resize(size, 0.)
  u.Resize(size);
  u.Clear();
  v.Resize(size);
  v.Clear();
  a.Resize(size);
  a.Clear();
  ut.Resize(size);
  ut.Clear();
  vt.Resize(size);
  vt.Clear();
  at.Resize(size);
  at.Clear();
  lambdaTrial = 1.0;
  lambdaConvg = 1.0;
/*
  //-------------------------------------------------------------------------
  // Apply initial conditions and gather displacements and accelerations
  //-------------------------------------------------------------------------
  theLoadCase->initialize();
  int i, n;
  n = pA->get_model()->get_model_nodes().size();
  for (i = 0; i < n; i++) {
    ModelNode* pModelNode = pA->get_model()->get_model_nodes()[i];
    IDContainer FTable = pModelNode->get_FTable();
    // Copy displacements
    pModelNode->zeroVector();
    pModelNode->add_uTrial();
    for (unsigned j = 0;j < FTable.size();j++)
      if ((FTable[j])>=0)
        u[FTable[j]]=pModelNode->get_vector()[j];
    // Copy velocities
    pModelNode->zeroVector();
    pModelNode->add_vTrial();
    for (unsigned j = 0;j < FTable.size();j++)
      if ((FTable[j])>=0)
        v[FTable[j]]=pModelNode->get_vector()[j];
  }
  //-------------------------------------------------------------------------
  // Apply loads
  //-------------------------------------------------------------------------
  pA->get_domain()->zeroLoads();
  theLoadCase->set_factor(1.0);
  theLoadCase->applyLoads(0.);

  //-------------------------------------------------------------------------
  // Build the system of equations M.a = R-K.u-C.v
  //-------------------------------------------------------------------------
  pA->get_soe()->zeroA();
  pA->get_soe()->zeroB();
  // Take contribution from Elements
  n = pA->get_model()->get_model_elements().size();
  for (i = 0; i < n; i++) {
    ModelElement* pModelElement = pA->get_model()->get_model_elements()[i];
    pModelElement->zeroMatrix();
    pModelElement->add_M();
    pA->get_soe()->insertMatrixIntoA(pModelElement->get_matrix(),
      pModelElement->get_FTable(), 1.0);

    pModelElement->zeroVector();
//    pModelElement->add_KuTrial(1.0);
//    pModelElement->add_CvTrial(1.0);
    pModelElement->add_Reff(1.0);
    pA->get_soe()->insertVectorIntoB(pModelElement->get_vector(),
      pModelElement->get_FTable(), -1.0);
  }
  // Take contribution from Nodes
  n = pA->get_model()->get_model_nodes().size();
  for (i = 0; i < n; i++) {
    ModelNode* p = pA->get_model()->get_model_nodes()[i];
    this->formNodalResidual(p);
    pA->get_soe()->insertVectorIntoB(p->get_vector(), p->get_FTable(), -1.0);
  }
  //-------------------------------------------------------------------------
  // Solve the system, update accelerations and commit
  //-------------------------------------------------------------------------
  pA->get_soe()->solve();
  a = pA->get_soe()->get_X();
  // Only accelerations are updated here, u and v are already updated from db
  pA->get_model()->set_trial_vecs(u, v, a);
  this->commit();
*/
}
/**
 * Forms the elemental tangent for a static analysis, element by element.
 * @param pModelElement A pointer to the ModelElement that is treated.
 * @return 0 if everything is ok.
 */
void TransientControl::FormElementalTangent(ModelElement* pModelElement)  {
  pModelElement->zeroMatrix();
  pModelElement->add_K(c[0]);
  pModelElement->add_M(c[1]);
  pModelElement->add_C(c[2]);
}
void TransientControl::FormElementalResidual(ModelElement* pModelElement,
                                             double /*time*/) {
  pModelElement->zeroVector();
  pModelElement->add_Reff(1.0);
}
void TransientControl::FormNodalResidual(ModelNode* pModelNode) {
  pModelNode->zeroVector();
  pModelNode->add_R(1.0);
}
void TransientControl::Commit() {
  pA->get_model()->Commit();
}
/**
 * Aborts step.
 */
void TransientControl::Rollback() {
  lambdaTrial = lambdaConvg;
}
void TransientControl::FormResidual(double fac)  {
  pA->get_soe()->zeroB();
  pA->get_domain()->zeroLoads();
  pA->get_domain()->ApplyLoads(fac, pA->get_domain()->get_time_curr());

  // Take contribution from Nodes
  for (unsigned i = 0; i < pA->get_model()->get_model_nodes().size(); i++) {
    ModelNode* p = pA->get_model()->get_model_nodes()[i];
    this->FormNodalResidual(p);
    pA->get_soe()->insertVectorIntoB(p->get_vector(), p->get_FTable(), -1.0);
  }
  // Take contribution from Elements
  for (unsigned i = 0; i < pA->get_model()->get_model_elements().size(); i++) {
    ModelElement* p = pA->get_model()->get_model_elements()[i];
    this->FormElementalResidual(p);
    pA->get_soe()->insertVectorIntoB(p->get_vector(), p->get_FTable(), -1.0);
  }
}
