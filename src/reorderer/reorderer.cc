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

#include "reorderer/reorderer.h"
#include <stdio.h>
#include "analysis/analysis.h"
#include "exception/sexception.h"
#include "model/model.h"
#include "model/model_element.h"
#include "model/model_node.h"

using std::cout;
using std::endl;

Reorderer::Reorderer() {
  pA->get_model()->set_reordered(false);
}
Reorderer::~Reorderer() {
}
int Reorderer::reorder() {
  if (pA->get_model()->isReordered()) return 0;
  std::vector<int> perm;
  // Get the permutation matrix
  if (this->get_perm(perm)>0) {
    // Reorder Model Nodes
    for (unsigned k = 0; k < pA->get_model()->get_model_nodes().size(); k++) {
      ModelNode* pModelNode = pA->get_model()->get_model_nodes()[k];
      for (unsigned i = 0; i < pModelNode->get_FTable().size(); i++) {
        if (pModelNode->get_FTable()[i] >= 0)
          pModelNode->set_FTable(i, perm[pModelNode->get_FTable()[i]]);
      }
    }
    // Reorder Model Elements
    for (unsigned k = 0; k < pA->get_model()->get_model_elements().size(); k++) {
      ModelElement* pModelElem = pA->get_model()->get_model_elements()[k];
      for (unsigned i = 0;i < pModelElem->get_FTable().size();i++) {
        if (pModelElem->get_FTable()[i] >= 0)
          pModelElem->set_FTable(i, perm[pModelElem->get_FTable()[i]]);
      }
    }
    printf("reo: Optimization returned successfully.\n");
  } else {
    printf("reo: Optimization failed.\n");
  }
  pA->get_model()->set_reordered(true);
  return 0;
}
