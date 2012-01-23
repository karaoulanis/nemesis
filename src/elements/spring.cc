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

#include "elements/spring.h"
#include <vector>
#include "group/group_data.h"
#include "main/nemesis_debug.h"
#include "material/spring_material.h"
#include "node/node.h"

/**
 * Default constructor.
 */
Spring::Spring()
  : dim_(0),
    mySpringMaterial(NULL) {
}

/**
 * Constructor.
 * Creates a Bar Element.
 */
Spring::Spring(int id,
               std::vector<Node*> nodes,
               SpringMaterial* material,
               int dim,
               double xp1,
               double xp2,
               double xp3,
               double yp1,
               double yp2,
               double yp3)
    : Element(id, nodes),
      dim_(dim) {
  // The dofs needed for this element
  myLocalNodalDofs.resize(dim_);
  for (int i = 0; i < dim_; i++) myLocalNodalDofs[i] = i;

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

  // Store material information
  myMaterial = material;
  mySpringMaterial = static_cast<SpringMaterial*>(myMaterial)->get_clone();

  // Find gap
  // gap = sqrt((x(1, 0)-x(0, 0))*(x(1, 0)-x(0, 0))
  //      +(x(1, 1)-x(0, 1))*(x(1, 1)-x(0, 1))
  //      +(x(1, 2)-x(0, 2))*(x(1, 2)-x(0, 2)));
  // Define transformations
  T.Resize(dim_, dim_, 0.);
  static Vector xp, yp, zp;
  switch (dim_) {
    case 1:
      T(0, 0)=1.;
      break;
    case 2:
      xp.Resize(2);
      yp.Resize(2);
      xp[0]= xp1;
      xp[1]= xp2;
      yp[0]=-xp2;
      yp[1]= xp1;
      xp.Normalize();
      yp.Normalize();
      T.AppendRow(xp, 0, 0);
      T.AppendRow(yp, 1, 0);
      break;
    case 3:
      xp.Resize(3);
      yp.Resize(3);
      zp.Resize(3);
      xp[0]=xp1;
      xp[1]=xp2;
      xp[2]=xp3;
      yp[0]=yp1;
      yp[1]=yp2;
      yp[2]=yp3;
      zp = Cross(xp, yp);
      yp = Cross(zp, xp);
      xp.Normalize();
      yp.Normalize();
      zp.Normalize();
      T.AppendRow(xp, 0, 0);
      T.AppendRow(yp, 1, 0);
      T.AppendRow(zp, 2, 0);
      break;
    default:
      break;
  }
}
/**
 * Destructor.
 */
Spring::~Spring() {
  delete mySpringMaterial;
}
void Spring::Update() {
  static Vector du(2*dim_);
  du = this->get_disp_incrm();
  report(du, "du");
  static Vector de(dim_, 0);
  de.Clear();
  for (int k = 0; k < dim_; k++)
    for (int i = 0; i < dim_; i++)
      de[i]+=T(i, k)*(du[k+dim_]-du[k]);
  // report(de, "de");
  mySpringMaterial->set_strain(de);
}
void Spring::Commit() {
  mySpringMaterial->Commit();
}


const Matrix& Spring::get_K() {
  Matrix& K=*myMatrix;
  K.Clear();
  const Matrix& Ct = mySpringMaterial->get_C();
  Matrix K0 = Transpose(T)*Ct*T;
  K.Append(K0, 0, 0, +1.);
  K.Append(K0, 0, 3, -1.);
  K.Append(K0, 3, 0, -1.);
  K.Append(K0, 3, 3, +1.);
  double facK = groupdata_->active ? groupdata_->factor_K: 1e-7;
  K*=facK;
  // report(K, "K", 14, 1);
  return K;

/*  Matrix& K=*myMatrix;
  K.Clear();
  Vector Ct = mySpringMaterial->get_C();
  for (int k = 0; k < nDim; k++) {
    for (int j = 0; j < nDim; j++) {
      for (int i = 0; i < nDim; i++) {
        K(i     , j     )=K(i     , j     )+T(k, i)*Ct[k]*T(k, j);
        K(i     , j+nDim)=K(i     , j+nDim)-T(k, i)*Ct[k]*T(k, j);
        K(i+nDim, j     )=K(i+nDim, j     )-T(k, i)*Ct[k]*T(k, j);
        K(i+nDim, j+nDim)=K(i+nDim, j+nDim)+T(k, i)*Ct[k]*T(k, j);
      }
    }
  }
  double facK = groupdata_->active_ ? groupdata_->factor_K_: 1e-7;
  K*=facK;
  return K;*/
}
const Matrix& Spring::get_M() {
  Matrix& M=*myMatrix;
  M.Clear();
  return M;
}
const Vector& Spring::get_R() {
  Vector& R=*myVector;
  R.Clear();
  // Quick return if inactive
  if (!(groupdata_->active)) {
    return R;
  }
  // Get factors
  double facS = groupdata_->factor_S;
  // double facG = groupdata_->factor_G;
  // double facP = groupdata_->factor_P;
  /// @todo check
  // double facG = myGroup->get_fac_G();
  // double facP = myGroup->get_fac_P();
  // R = Fint - Fext
  Vector sigma = mySpringMaterial->get_stress();
  for (int k = 0; k < dim_; k++) {
    for (int i = 0; i < dim_; i++) {
      double d = facS*T(k, i)*sigma[k];
      R[i     ]+=-d;
      R[i+dim_]+= d;
    }
  }
  report(R, "R");
  return R;
}
const Vector& Spring::get_Reff() {
  /// @todo: form Reff
  Vector& Reff=*myVector;
  Reff.Clear();
  return Reff;
}
const Vector& Spring::get_Rgrad() {
  /// @todo: form Rgrad
  Vector& Rgrad=*myVector;
  Rgrad.Clear();
  return Rgrad;
}
void Spring::recoverStresses() {
  /// @todo Stresses from bar to nodes
  static Vector s(6);
  s.Clear();
  s.Append(mySpringMaterial->get_stress(), 0);
  nodes_[0]->addStress(s);
  nodes_[1]->addStress(s);
}
