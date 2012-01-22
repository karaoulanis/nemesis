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

#include "material/creep.h"
#include <sstream>
#include <string>

Matrix Creep::C(6, 6, 0.);

Creep::Creep()
    : myElastic(0),
      A(0.),
      n(0.),
      k(0.) {
}

Creep::Creep(int id, MultiaxialMaterial* elastic, double A, double n, double k)
    : MultiaxialMaterial(id, 0., 0.),
      A(0.),
      n(0.),
      k(0.) {
  // Get the elastic part
  myElastic = elastic->get_clone();
  MatParams[ 0]=A;
  MatParams[ 1]=n;
  MatParams[ 2]=k;
  MatParams[30]=myElastic->get_param(30);
  MatParams[31]=myElastic->get_param(31);

  // Material state
  eCTrial.Resize(6, 0.);
  eCConvg.Resize(6, 0.);
}
Creep::~Creep() {
  delete myElastic;
}

/**
 * Update stresses given a total strain increment.
 * @param De Vector containing total strain increment.
 */
void Creep::set_strain(const Vector& /*De*/) {
  /// @todo
  // double Dt = pD->get_time_incr();
  // eTrial = eTotal+De;
  // // sTrial = sConvg+(myElastic->get_C())*De;
  // sTrial = (myElastic->get_C())*(eTotal-eCConvg);
  // double d = eCConvg.twonorm()*sqrt(3./2.);
  // Vector beta(6, 0);
  // beta = k*pow(A, 1./k)*d*sTrial;
  // eCTrial = eCConvg+Dt*beta;
  // sTrial=(myElastic->get_C())*(eTotal-eCTrial-Dt*beta);
}

/**
 * Commit material state.
 */
void Creep::Commit() {
  eTotal = eTrial;  /// @todo
  eCConvg = eCTrial;
  sConvg = sTrial;
}

/**
 * Get tangent material matrix.
 * @todo fill it
 * @return A reference to the tangent material matrix.
 */
const Matrix& Creep::get_C() {
  return myElastic->get_C();
}

MultiaxialMaterial* Creep::get_clone() {
  // Material parameters
  int myID    = this->get_id();
  double A    = MatParams[ 0];
  double n    = MatParams[ 1];
  double k    = MatParams[ 2];
  // Create clone and return
  Creep* newClone = new Creep(myID, myElastic, A, n, k);
  return newClone;
}

void Creep::Save(std::ostream* os) {
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
