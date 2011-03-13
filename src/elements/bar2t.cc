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

#include "elements/bar2t.h"
#include "loadcase/group_data.h"
#include "material/uniaxial_material.h"

/**
 * Default constructor.
 */
Bar2t::Bar2t()   {
}
/**
 * Constructor.
 * Creates a Bar2t Element. At this point only the id's of the nodes and the
 * Material are stored. Only when this element is added to the Domain, the
 * pointers concerning the Nodes and the Material are initiated.
 */
Bar2t::Bar2t(int ID, int Node_1, int Node_2, int matID, int iSecID, int jSecID)
  :Bar(ID, Node_1, Node_2, matID, iSecID, jSecID) {
  myTag = TAG_ELEM_BAR_2D_TOTAL_LAGRANGIAN;
  cosX.resize(nDim);
  cosX.clear();
}
/**
 * Destructor.
 */
Bar2t::~Bar2t() {
}
const Matrix& Bar2t::get_K() {
  Matrix& K=*myMatrix;

  // Directional cosines, L^2 and strain (in the current configuration)
  static Vector du(2*nDim);
  du = this->get_disp_trial();
  for (int i = 0;i < nDim;i++) cosX[i]=(x(1, i)-x(0, i)+du[i+nDim]-du[i])/L0;
  double L2 = 0;
  for (int i = 0;i < nDim;i++)
    L2+=(x(1, i)-x(0, i)+du[i+nDim]-du[i])*(x(1, i)-x(0, i)+du[i+nDim]-du[i]);
  double e=(L2-L0*L0)/(2*L0*L0);

  // Matrix KM
  K.clear();
  double E = myUniMaterial->get_C();
  double K0 = E*A0/L0;
  double d;
  for (int i = 0;i < nDim;i++)
    for (int j = 0; j < nDim; j++) {
      d = cosX[i]*cosX[j]*K0;
      K(i     , j)      = d;
      K(i+nDim, j)      =-d;
      K(i     , j+nDim) =-d;
      K(i+nDim, j+nDim) = d;
    }
  // Matrix KG
  K0 = E*A0*e/L0;
  int ii;
  for (int i = 0; i < nDim; i++) {
    i < nDim ? ii = i+nDim : ii = i-nDim;
    K(i, i) += K0;
    K(i, ii)+=-K0;
  }

  /// @todo Add nodal transformations
  double facK = groupdata_->active_ ? groupdata_->factor_K_: 1e-7;
  K*=facK;
  return K;
}
const Vector& Bar2t::get_R() {
  Vector& R=*myVector;
  R.clear();

  // Factors
  if (!(groupdata_->active_))  return R;
  double facS = groupdata_->factor_S_;
  double facG = groupdata_->factor_G_;
  double facP = groupdata_->factor_P_;

  // Directional cosines, L^2 and strain (in the current configuration)
  static Vector du(2*nDim);
  du = this->get_disp_trial();
  for (int i = 0;i < nDim;i++) cosX[i]=(x(1, i)-x(0, i)+du[i+nDim]-du[i])/L0;
  double L2 = 0;
  for (int i = 0;i < nDim;i++)
    L2+=(x(1, i)-x(0, i)+du[i+nDim]-du[i])*(x(1, i)-x(0, i)+du[i+nDim]-du[i]);
  double e=(L2-L0*L0)/(2*L0*L0);
  double E = myUniMaterial->get_C();
  double N = E*A0*e;
  for (int i = 0; i < nDim; i++) {
    double d = facS*cosX[i]*N;
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
