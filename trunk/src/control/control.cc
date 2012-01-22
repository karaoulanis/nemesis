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

#include "control/control.h"
#include "analysis/analysis.h"
#include "model/model.h"
#include "model/model_element.h"
#include "model/model_node.h"
#include "soe/soe.h"

Control::Control()
    : lambdaTrial(0.0),
      lambdaConvg(0.0),
      DLambda(0.0),
      dLambda(0.0),
      qRef(0) {
}

Control::~Control() {
}

void Control::FormTangent() {
  pA->get_soe()->zeroA();
  int n = pA->get_model()->get_model_elements().size();
  for (int i = 0; i < n; i++) {
    ModelElement* p = pA->get_model()->get_model_elements()[i];
    this->FormElementalTangent(p);
    pA->get_soe()->InsertMatrixIntoA(p->get_matrix(), p->get_FTable(), 1.0);
  }
}

double Control::get_lambda() {
  return lambdaTrial;
}
