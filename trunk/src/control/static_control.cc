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

#include "control/static_control.h"
#include "analysis/analysis.h"
#include "domain/domain.h"
#include "model/model_element.h"
#include "model/model_node.h"
#include "soe/soe.h"

/**
 * Constructor.
 */
StaticControl::StaticControl()
    : Delta0(0.),
      minDelta(0.),
      maxDelta(fabs(0.)),
      du(0),
      duT(0),
      duBar(0),
      Du(0),
      Io(0),
      Id(0),
      nExp(0.),
      Dt(0.) {
}
/**
 * Constructor.
 * \param D0      Initial Delta.
 * \param minD      Lower bound for Delta.
 * \param maxD      Upper bound for Delta.
 * \param n       Exponent parameter.
 * \param IterDesired Desired number of iterations.
 * \param DeltaTime     Timestep for viscoplastic solutions.
 */
StaticControl::StaticControl(double D0, double minD, double maxD,
                   int IterDesired, double n, double DeltaTime)
    : Delta0(D0),
      minDelta(fabs(minD)),
      maxDelta(fabs(maxD)),
      du(0),
      duT(0),
      duBar(0),
      Du(0),
      Io(IterDesired),
      Id(IterDesired),
      nExp(n),
      Dt(DeltaTime) {
  DLambda = D0;
  if (minDelta>maxDelta)
    throw SException("[nemesis:%d] %s", 9999, "'max' is less than 'min'.");
  if (IterDesired <= 0)
    throw SException("[nemesis:%d] %s", 9999, "'Id' should be greater than 0.");
}
/**
 * Destructor.
 */
StaticControl::~StaticControl() {
}
/**
 * Initializes the vectors used in the Static Control. Also finds the
 * vector of external forces and initializes the norm.
 * @return 0 if everything is ok.
 */
void StaticControl::Init() {
  // Check the size
  int size = pA->get_model()->get_num_eqns();
  /// @todo resize(size, 0.)
  qRef.Resize(size);
  qRef.Clear();
  Du.Resize(size);
  Du.Clear();
  du.Resize(size);
  du.Clear();
  duBar.Resize(size);
  duBar.Clear();
  duT.Resize(size);
  duT.Clear();

  this->FormResidual(1.0);
  qRef = pA->get_soe()->get_B();

  lambdaTrial = 0;
  lambdaConvg = 0;
  dLambda = 0;
  pA->get_domain()->incTime(Dt);
}
/**
 * Forms the elemental tangent for a static analysis, element by element.
 * @param pModelElement A pointer to the ModelElement that is treated.
 * @return 0 if everything is ok.
 */
void StaticControl::FormElementalTangent(ModelElement* pModelElement)  {
  pModelElement->zeroMatrix();
  pModelElement->add_K(1.0);
}
/**
 * Forms the elemental residual for a static analysis, element by element.
 * @param pModelElement A pointer to the ModelElement that is treated.
 * @param time Current(?) time.
 * @return 0 if everything is ok.
 */
void StaticControl::FormElementalResidual(ModelElement* pModelElement,
                                          double /*time*/) {
  pModelElement->zeroVector();
  pModelElement->add_R(1.0);
}
/**
 * Forms the nodal tangent for a static analysis, node by node.
 * @param pModelNode A pointer to the ModelNode that is treated.
 * @return 0 if everything is ok.
 */
void StaticControl::FormNodalResidual(ModelNode* pModelNode) {
  pModelNode->zeroVector();
  pModelNode->add_R(1.0);
}
void StaticControl::FormResidual(double fac)   {
  pA->get_soe()->zeroB();
  pA->get_domain()->zeroLoads();
  pA->get_domain()->ApplyLoads(fac, 0.);

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
/**
 * Commits displacements.
 * After a successfull step the displacements in the domain are set to be
 * converged. Time and lambda are committed to Domain.
 * @return 0 if everything is ok.
 */
void StaticControl::Commit() {
  lambdaConvg = lambdaTrial;
  pA->get_domain()->set_lambda(lambdaConvg);
  pA->get_domain()->Commit(); /// @todo this commits only domains time!
  pA->get_model()->Commit();
}
/**
 * Aborts step.
 * After a unsuccessfull step the displacement in the domain are set to the
 * previous converged step.
 * @return 0 if everything is ok.
 */
void StaticControl::Rollback() {
  // pA->get_model()->rollBack();
  lambdaTrial = lambdaConvg;
}

