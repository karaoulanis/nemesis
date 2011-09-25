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

#include "elements/bar.h"
#include <vector>
#include "crosssection/cross_section.h"
#include "material/uniaxial_material.h"
#include "node/node.h"

/**
 * Default constructor.
 */
Bar::Bar()
    : dim_(0),
      iSection(0),
      jSection(0),
      L0(0.),
      A0(0.),
      cosX(0),
      myUniMaterial(0) {
}
/**
 * Constructor.
 * Creates a Bar Element.
 */
Bar::Bar(int id, std::vector<Node*> nodes, UniaxialMaterial* material,
         CrossSection* iSec, CrossSection* jSec, int dim)
    : Element(id, nodes),
      dim_(dim),
      iSection(iSec),
      jSection(jSec) {

  // Store material information
  myUniMaterial = material->get_clone();

  // Get nodal data
  myNodalIDs.resize(2);
  myNodalIDs[0] = nodes_[0]->get_id();
  myNodalIDs[1] = nodes_[1]->get_id();

  // Set local nodal dofs
  myLocalNodalDofs.resize(dim_);
  for (int i = 0; i < dim_; i++) {
    myLocalNodalDofs[i] = i;
  }

  // Handle common info: Start -------------------------------------------------
  // Find own matrix and vector
  myMatrix = theStaticMatrices[2*dim_];
  myVector = theStaticVectors[2*dim_];
  // Get nodal coordinates
  /// @todo replace with const references
  x.Resize(2, 3);
  for (int i = 0; i < 2; i++) {
    x(i, 0) = nodes_[i]->get_x1();
    x(i, 1) = nodes_[i]->get_x2();
    x(i, 2) = nodes_[i]->get_x3();
  }
  // Inform the nodes that the corresponding Dof's must be activated
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < dim_; j++)
      nodes_[i]->addDofToNode(myLocalNodalDofs[j]);
  // Load vector
  P.Resize(2*dim_, 0.);
  // Self weight
  G.Resize(2*dim_, 0.);
  this->AssignGravityLoads();
  // Handle common info: End ---------------------------------------------------

  // Find length
  L0 = 0.;
  for (int i = 0; i < dim_; i++) {
    L0 += (x(1, i)-x(0, i))*(x(1, i)-x(0, i));
  }
  L0 = sqrt(L0);
  if (num::tiny(L0)) {
    throw SException("[nemesis:%d] %s", 9999, "Zero length bar is not allowed");
  }
  // Retrieve the CrossSection pointers and get A0
  A0 = 0.5*(iSection->get_A()+jSection->get_A());
}


Bar::~Bar() {
  delete myUniMaterial;
}


void Bar::update() {
  static Vector du(2*dim_);
  du = this->get_disp_incrm();
    double dL = 0;
  for (int i = 0;i < dim_; i++) dL+=(du[i+dim_]-du[i])*cosX[i];
    double de = dL/L0;
  myUniMaterial->set_strain(de);
}


void Bar::commit() {
  myUniMaterial->commit();
}


const Matrix& Bar::get_M() {
  Matrix& M=*myMatrix;
  M.Clear();
  double rho = myUniMaterial->get_rho();
  double mass = 0.5*L0*A0*rho;
  for (int i = 0; i < dim_; i++) {
    M(i     , i)      = mass;
    M(i+dim_, i+dim_) = mass;
  }
  return M;
}


const Vector& Bar::get_Reff() {
  /// @todo: problem with memory sharing otherwise
  Vector velc = this->get_velc_trial();
  Vector& Reff=*myVector;
  // +Fint-Fext
  this->get_R();
  // +M*aTrial
  double rho = myUniMaterial->get_rho();
  double mass = 0.5*rho*L0;
  const Vector& a0 = nodes_[0]->get_accl_trial();
  const Vector& a1 = nodes_[1]->get_accl_trial();
  for (int i = 0;i < dim_;i++) {
      Reff[i]      += mass*a0[i];
      Reff[i+dim_] += mass*a1[i];
  }
  Reff+=(this->get_C())*(velc);
  return Reff;
}


void Bar::recoverStresses() {
  /// @todo Stresses from bar to nodes
  static Vector s(6);
  s[0]=myUniMaterial->get_stress();
  nodes_[0]->addStress(s);
  nodes_[1]->addStress(s);
}
