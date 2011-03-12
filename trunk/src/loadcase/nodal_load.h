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

#ifndef SRC_LOADCASE_NODAL_LOAD_H_
#define SRC_LOADCASE_NODAL_LOAD_H_

#include "loadcase/load.h"

class Node;

class NodalLoad: public Load {
 public:
  /**
   * Default constructor.
   */
  NodalLoad();

  /**
   * Constructor.
   * @param node The node where the Load is apllied.
   * @param dof The corresponding dof (ux = 1, uy = 2, ...)
   */
  NodalLoad(Node* node, int dof);

  /**
   * Virtual destructor.
   */
  virtual ~NodalLoad();

  /**
   * Apply load to its Node.
   * @param factor Factor multiplied by the Load value.
   * @param time The time at which Load value is evaluated.
   */
  void Apply(double factor, double time);

  /**
   * Find Load value at given time.
   * @param time The time at which Load value is evaluated.
   */
  virtual double GetValue(double time)=0;

 protected:
  Node* node_;
  int dof_;
};

#endif  // SRC_LOADCASE_NODAL_LOAD_H_
