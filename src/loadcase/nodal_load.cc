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

#include "loadcase/nodal_load.h"
#include "node/node.h"

NodalLoad::NodalLoad() {
}

NodalLoad::~NodalLoad() {
}

NodalLoad::NodalLoad(Node* node, int dof)
:Load() {
  node_ = node;
  dof_  = dof - 1;
  ///@todo replace NodalLoad::get_activated_dof(dof_) by IsDofActive()
  ///@todo Give SException the right id and not 9999.
  if (node_->get_activated_dof(dof_) < 0)
    throw SException("[nemesis:%d] %s", 9999, "Dof is not activated.");
}

void NodalLoad::Apply(double factor, double time) {
  node_->addLoad(dof_, this->GetValue(time), factor);
}
