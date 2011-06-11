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

#include "elements/bar2s.h"
#include <vector>
#include "group/group_data.h"
#include "material/uniaxial_material.h"

/**
 * Default constructor.
 */
Bar2s::Bar2s()   {
}

/**
 * Constructor.
 * Creates a Bar2s Element.
 */
Bar2s::Bar2s(int id, std::vector<Node*> nodes, UniaxialMaterial* material,
         CrossSection* iSec, CrossSection* jSec, int dim)
    : Bar(id, nodes, material, iSec, jSec, dim) {
  // Find directional cosines
  cosX.resize(dim_);
  cosX.clear();
  for (int i = 0;i < dim_; i++) cosX[i]=(x(1, i)-x(0, i))/L0;
}
/**
 * Destructor.
 */
Bar2s::~Bar2s()  {
}
const Matrix& Bar2s::get_K() {
  Matrix& K=*myMatrix;
  double E = myUniMaterial->get_C();
  double K0 = E*A0/L0;
  double d;
  for (int i = 0;i < dim_; i++)
    for (int j = 0; j < dim_; j++) {
      d = cosX[i]*cosX[j]*K0;
      K(i     , j)      =  d;
      K(i+dim_, j)      = -d;
      K(i     , j+dim_) = -d;
      K(i+dim_, j+dim_) =  d;
    }
  /// @todo Add nodal transformations
  double facK = groupdata_->active ? groupdata_->factor_K : 1e-7;
  K*=facK;
  return K;
}
const Vector& Bar2s::get_Rgrad() {
  double E = myUniMaterial->get_C();
  Matrix& K=*myMatrix;
  double K0;
  switch (activeParameter) {
  case 1:  K0 = A0/L0; break;  // E is the active parameter
  case 2:  K0= E/L0; break;  // A is the active parameter
  case 0:
  default: K0 = 0.;
  }
  double d;
  for (int i = 0;i < dim_;i++)
    for (int j = 0; j < dim_; j++) {
      d = cosX[i]*cosX[j]*K0;
      K(i     , j)      = d;
      K(i+dim_, j)      =-d;
      K(i     , j+dim_) =-d;
      K(i+dim_, j+dim_) = d;
    }
  /// @todo Add nodal transformations
  double facK = groupdata_->active ? groupdata_->factor_K : 1e-7;
  K*=facK;
  *myVector = K*(this->get_disp_convg());
  return *myVector;
}
const Vector& Bar2s::get_R() {
  Vector& R=*myVector;
  R.clear();

  // Quick return if element not active
  if (!(groupdata_->active)) {
    return R;
  }

  // Factors
  double facS = groupdata_->factor_S;
  double facG = groupdata_->factor_G;
  double facP = groupdata_->factor_P;

  // R = Fint - Fext
  double F0 = A0*(myUniMaterial->get_stress());
  for (int i = 0; i < dim_; i++) {
    double d = facS*cosX[i]*F0;
    R[i]      =-d - facP*P[i];
    R[i+dim_] = d - facP*P[i+dim_];
  }
  // Self-weigth (only in y todo: check this)
  if (!num::tiny(facG) && dim_>1) {
    double b=-facG*0.5*(myUniMaterial->get_rho())*A0*L0;
    R[1]-=b;
    R[1+dim_]-=b;
  }
  /// @todo Add nodal transformations
  return R;
}
