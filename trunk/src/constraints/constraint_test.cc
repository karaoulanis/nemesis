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

// Included files
#include "constraints/constraint.h"
#include <gtest/gtest.h>
#include <exception>
#include <domain/domain.h>
#include <node/node.h>

class ConstraintTest : public ::testing::Test {
 protected:
  Domain d;
  Node* n;
  Constraint* c;
  CDof cdof;
  
  virtual void SetUp() {
    d.set_dim(2);
    n = new Node(1, 0., 0.);
    d.add(d.get_nodes(), n);
    c = new Constraint();
    // c->set_val(2.);
  }

  virtual void TearDown() {
    // delete n;
    delete c;
  }
};

TEST_F(ConstraintTest, DofNotActivated) {
  EXPECT_THROW({
    // Take the already existing node 
    Node* n = d.get<Node>(d.get_nodes(), 1);
    // Try to add a constraint to a not activated dof
    c->set_cdof(n, 1, 8.);
  }, SException);
}
