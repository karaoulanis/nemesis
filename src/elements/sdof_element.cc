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

#include "elements/sdof_element.h"
#include <vector>
#include "group/group_data.h"
#include "node/node.h"

/**
 * Default constructor.
 */
SDofElement::SDofElement()
    : mySDofMaterial(0) {
}

/**
 * Constructor.
 */
SDofElement::SDofElement(int id, std::vector<Node*> nodes, int dof,
                         SDofMaterial* material)
    : Element(id, nodes),
      mySDofMaterial(material) {
  // The dofs needed for this element
  myNodalIDs.resize(1);
  myNodalIDs[0] = nodes_[0]->get_id();
  // Local dofs
  myLocalNodalDofs.resize(1);
  myLocalNodalDofs[0] = dof-1;
  // Handle common info: Start -------------------------------------------------
  // Find own matrix and vector
  myMatrix = theStaticMatrices[myLocalNodalDofs.size()*nodes_.size()];
  myVector = theStaticVectors[myLocalNodalDofs.size()*nodes_.size()];
  // Get nodal coordinates
  /// @todo replace with const references
  x.Resize(nodes_.size(), myLocalNodalDofs.size());
  for (unsigned i = 0; i < nodes_.size(); i++) {
    x(i, 0) = nodes_[i]->get_x1();
    x(i, 1) = nodes_[i]->get_x2();
    x(i, 2) = nodes_[i]->get_x3();
  }
  // Inform the nodes that the corresponding Dof's must be activated
  for (unsigned i = 0; i < nodes_.size(); i++) {
    for (unsigned j = 0; j < myLocalNodalDofs.size() ; j++) {
      nodes_[i]->addDofToNode(myLocalNodalDofs[j]);
    }
  }
  // Load vector
  P.resize(myLocalNodalDofs.size()*nodes_.size(), 0.);
  // Self weight
  G.resize(myLocalNodalDofs.size()*nodes_.size(), 0.);
  this->AssignGravityLoads();
  // Handle common info: End ---------------------------------------------------
}

SDofElement::~SDofElement() {
}

const Matrix& SDofElement::get_K() {
  Matrix& K=*myMatrix;
  double facK = groupdata_->active ? groupdata_->factor_K : 1e-7;
  K(0, 0)=facK*(mySDofMaterial->get_param(0));
  return K;
}

const Matrix& SDofElement::get_M() {
  Matrix& M=*myMatrix;
  M(0, 0)=mySDofMaterial->get_rho();
  return M;
}

const Vector& SDofElement::get_R() {
  Vector& R=*myVector;
  R.clear();
  return R;
}

bool SDofElement::checkIfAllows(FEObject* /*f*/) {
  return true;
}
