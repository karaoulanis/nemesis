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

#include "group/group.h"

Group::Group() {
}
Group::Group(int ID)
:DomainObject(ID)
{ 
  active = true;
  facK = 1.;
  facS = 1.;
  facG = 1.;
  facP = 1.;
}
Group::~Group() {
}
void Group::set_default() {
  facK = 1.;
  facS = 1.;
  facG = 1.;
  facP = 1.;
}
void Group::set_state(GroupState* g) {
  active = g->get_active();
  facK = g->get_fac_K();
  facS = g->get_fac_S();
  facG = g->get_fac_G();
  facP = g->get_fac_P();
}
