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

#include "material/spring_elastic.h"

SpringElastic::SpringElastic() {
}
SpringElastic::SpringElastic(int ID, double Kn, double Ks2, double Ks3)
:SpringMaterial(ID) {
  // Material parameters
  MatParams[0]=Kn;
  MatParams[1]=Ks2;
  MatParams[2]=Ks3;
  myTag = TAG_MATERIAL_SPRING;

  // Initialize Ct
  Ct.clear();
  Ct(0, 0)=Kn;
  Ct(1, 1)=Ks2;
  Ct(2, 2)=Ks3;
}
SpringMaterial* SpringElastic::get_clone() {
  // Material parameters
  double Kn  =MatParams[0];
  double Ks2 =MatParams[1];
  double Ks3 =MatParams[2];
  // Create clone and return
  SpringMaterial* clone = new SpringElastic(myID, Kn, Ks2, Ks3);
  return clone;
}
void SpringElastic::set_strain(const Vector& De) {
  double Kn  =MatParams[0];
  double Ks2 =MatParams[1];
  double Ks3 =MatParams[2];
  eTrial = eTotal+De;
  switch (nDim) {
  case 3:
    sTrial[2]=Ks3*eTrial[2];
  case 2:
    sTrial[1]=Ks2*eTrial[1];
  case 1:
    sTrial[0]=Kn*eTrial[0];
  default:
    break;
  }
}
