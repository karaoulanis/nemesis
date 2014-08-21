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

#include "material/multiaxial_elastic.h"
#include <sstream>
#include <string>

MultiaxialElastic::MultiaxialElastic() {
}
MultiaxialElastic::MultiaxialElastic(int ID, double E, double nu, double rho,
                                     double aT, double kx, double ky, double kz)
:MultiaxialMaterial(ID, rho, aT) {
  // Material parameters
  MatParams[0]=E;
  MatParams[1]=nu;
  MatParams[2]=kx;
  MatParams[3]=ky;
  MatParams[4]=kz;
}
MultiaxialMaterial* MultiaxialElastic::get_clone() {
  // Material parameters
  double E   =MatParams[ 0];
  double nu  =MatParams[ 1];
  double kx  =MatParams[ 2];
  double ky  =MatParams[ 3];
  double kz  =MatParams[ 4];
  double rho =MatParams[30];
  double aT  =MatParams[31];
  // Create clone and return
  MultiaxialElastic* clone =
    new MultiaxialElastic(id_, E, nu, rho, aT, kx, ky, kz);
  return clone;
}
void MultiaxialElastic::set_strain(const Vector& De, const double Dt) {
  eTrial = eTotal+De;
  // sTrial=(this->get_C())*eTrial;
  sTrial = sConvg+(this->get_C())*De;
}
const Matrix& MultiaxialElastic::get_C() {
  C.Clear();
  // Material parameters
  double kx  = MatParams[ 2];
  double ky  = MatParams[ 3];
  double kz  = MatParams[ 4];
  double E   = MatParams[ 0]+(kx*x+ky*y+kz*z);
  double nu  = MatParams[ 1];
  // Find and return C
  double Em = E/((1.+nu)*(1.-2*nu));
  C(0, 0) = Em*(1.-nu);
  C(0, 1) = Em*nu;
  C(0, 2) = Em*nu;
  C(1, 0) = Em*nu;
  C(1, 1) = Em*(1.-nu);
  C(1, 2) = Em*nu;
  C(2, 0) = Em*nu;
  C(2, 1) = Em*nu;
  C(2, 2) = Em*(1.-nu);
  C(3, 3) = Em*0.5*(1.-2*nu);
  C(4, 4) = Em*0.5*(1.-2*nu);
  C(5, 5) = Em*0.5*(1.-2*nu);
  return C;
}
void MultiaxialElastic::Commit() {
  eTotal = eTrial;
  sConvg = sTrial;
}

void MultiaxialElastic::Save(std::ostream* os) {
  // start saving
  (*os) << "{";
  (*os) << "\"data\":{";
  (*os) << "\"sigm\":"    << sConvg << ',';
  (*os) << "\"epst\":"    << eTotal << ',';
//  (*os) << "\"epsp\":"    << ePConvg << ',';
  (*os) << "\"epsv\":"    << eTotal[0]+eTotal[1]+eTotal[2] << ',';
  (*os) << "\"p\":"       << sConvg.p() << ',';
  (*os) << "\"q\":"       << sConvg.q();
  (*os) << "}";
  // finalize
  (*os) << "}";
}
