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

#include "material/creep.h"

Matrix Creep::C(6, 6, 0.);

Creep::Creep() {
}
Creep::Creep(int ID, int elasticID, double A, double n, double k)
:MultiaxialMaterial(ID, 0., 0.) {
  // Get the elastic part
  Material* p = pD->get < Material>(pD->get_materials(), elasticID);
  if (p->get_tag() != TAG_MATERIAL_MULTIAXIAL_ELASTIC)
    throw SException("[nemesis:%d] %s", 9999,
                     "Multiaxial elastic material expected.");
  myElastic = static_cast < MultiaxialMaterial*>(p)->get_clone();
  MatParams[ 0]=A;
  MatParams[ 1]=n;
  MatParams[ 2]=k;
  MatParams[30]=myElastic->get_param(30);
  MatParams[31]=myElastic->get_param(31);

  // Material state
  eCTrial.resize(6, 0.);
  eCConvg.resize(6, 0.);
}
Creep::~Creep() {
  delete myElastic;
}
/**
 * Update stresses given a total strain increment.
 * @param De Vector containing total strain increment.
 */ 
void Creep::set_strain(const Vector& De) {
  double Dt = pD->get_time_incr();
  eTrial = eTotal+De;
  // sTrial = sConvg+(myElastic->get_C())*De;
  sTrial = (myElastic->get_C())*(eTotal-eCConvg);
  double d = eCConvg.twonorm()*sqrt(3./2.);
  Vector beta(6, 0);
  beta = k*pow(A, 1./k)*d*sTrial;
  eCTrial = eCConvg+Dt*beta;
  sTrial=(myElastic->get_C())*(eTotal-eCTrial-Dt*beta);
}
/**
 * Commit material state.
 */
void Creep::commit() {
  eTotal = eTrial;  /// @todo
  eCConvg = eCTrial;
  sConvg = sTrial;
  this->track();
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
  int elID    = myElastic->get_id();
  double A    = MatParams[ 0];
  double n    = MatParams[ 1];
  double k    = MatParams[ 2];
  // Create clone and return
  Creep* newClone = new Creep(myID, elID, A, n, k);
  return newClone;
}

/**
 * Add a record to the tracker.
 * If \a myTracker pointer is null (no tracker is added) just return.
 * Otherwise gather info and send them to the tracker.
 * The domain should be already updated!
 */
void Creep::track() {
  if (myTracker == 0) return;
  ostringstream s;
  s << "DATA "  << ' ';
  s << "sigm "  << ' ' <<sConvg;
  s << "epst "  << ' ' <<eTotal;
  s << "epsp "  << ' ' <<eCConvg;
  s << "epsv "  << 1020 << ' ' <<eTotal[0]+eTotal[1]+eTotal[2] << ' ';
  s << "p "     << 1020 << ' ' <<sConvg.p() << ' ';
  s << "q "     << 1020 << ' ' <<sConvg.q() << ' ';
  s << "END " <<' ';
  myTracker->track(pD->get_lambda(), pD->get_time_curr(), s.str());
}
