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

#include "material/plane_stress.h"

PlaneStress::PlaneStress() {
}
///@todo 0.
PlaneStress::PlaneStress(int ID, double E, double nu, double rho, double aT)
:MultiaxialElastic(ID, E, nu, rho, aT, 0., 0., 0.) {
}
MultiaxialMaterial* PlaneStress::getClone() {
  // Material parameters
  double E   =MatParams[ 0];
  double nu  =MatParams[ 1];
  double rho =MatParams[30];
  double aT  =MatParams[31];
  // Create clone and return
  PlaneStress* newClone = new PlaneStress(myID, E, nu, rho, aT);
  return newClone;
}
const Matrix& PlaneStress::getC() {
  C.clear();
  // Material parameters
  double E   =MatParams[ 0];
  double nu  =MatParams[ 1];
  // Find and return C
  double Em = E/(1-nu*nu);
  C(0, 0)=Em;
  C(0, 1)=Em*nu;
  C(1, 0)=Em*nu;
  C(1, 1)=Em;
  C(3, 3)=Em*0.5*(1.-nu);
  return C;
}
