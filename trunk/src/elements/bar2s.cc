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

#include "elements/bar2s.h"
#include "loadcase/group_data.h"
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
Bar2s::Bar2s(int ID, int Node_1, int Node_2, int matID,
             CrossSection* iSec, CrossSection* jSec)
  :Bar(ID, Node_1, Node_2, matID, iSec, jSec) {
  myTag = TAG_ELEM_BAR_2D_GEOMETRICALLY_LINEAR;
  // Find directional cosines
  cosX.resize(nDim);
  cosX.clear();
  for (int i = 0;i < nDim;i++) cosX[i]=(x(1, i)-x(0, i))/L0;
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
  for (int i = 0;i < nDim;i++)
    for (int j = 0; j < nDim; j++) {
      d = cosX[i]*cosX[j]*K0;
      K(i     , j)      =  d;
      K(i+nDim, j)      = -d;
      K(i     , j+nDim) = -d;
      K(i+nDim, j+nDim) =  d;
    }
  /// @todo Add nodal transformations
  double facK = groupdata_->active_ ? groupdata_->factor_K_: 1e-7;
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
  for (int i = 0;i < nDim;i++)
    for (int j = 0; j < nDim; j++) {
      d = cosX[i]*cosX[j]*K0;
      K(i     , j)      = d;
      K(i+nDim, j)      =-d;
      K(i     , j+nDim) =-d;
      K(i+nDim, j+nDim) = d;
    }
  /// @todo Add nodal transformations
  double facK = groupdata_->active_ ? groupdata_->factor_K_: 1e-7;
  K*=facK;
  *myVector = K*(this->get_disp_convg());
  return *myVector;
}
const Vector& Bar2s::get_R() {
  Vector& R=*myVector;
  R.clear();

  // Quick return if element not active
  if (!(groupdata_->active_))  return R;

  // Factors
  double facS = groupdata_->factor_S_;
  double facG = groupdata_->factor_G_;
  double facP = groupdata_->factor_P_;

  // R = Fint - Fext
  double F0 = A0*(myUniMaterial->get_stress());
  for (int i = 0; i < nDim; i++) {
    double d = facS*cosX[i]*F0;
    R[i]      =-d - facP*P[i];
    R[i+nDim] = d - facP*P[i+nDim];
  }
  // Self-weigth (only in y todo: check this)
  if (!num::tiny(facG) && nDim>1) {
    double b=-facG*0.5*(myUniMaterial->get_rho())*A0*L0;
    R[1]-=b;
    R[1+nDim]-=b;
  }
  /// @todo Add nodal transformations
  return R;
}
