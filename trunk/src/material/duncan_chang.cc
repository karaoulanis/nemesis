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

#include "material/duncan_chang.h"

DuncanChang::DuncanChang() {
}
DuncanChang::DuncanChang(int ID, double E, double nu, double c, double phi,
             double m, double Rf, double pa, double rho, double aT)
:MultiaxialMaterial(ID, rho, aT) {
  // Material parameters
  MatParams[0]=E;
  MatParams[1]=nu;
  MatParams[2]=m;
  MatParams[3]=c;
  MatParams[4]=phi;
  MatParams[5]=Rf;
  MatParams[6]=pa;
  // Material tag
  myTag = TAG_MATERIAL_MULTIAXIAL_ELASTIC;
}
MultiaxialMaterial* DuncanChang::get_clone() {
  // Material parameters
  double E   = MatParams[ 0];
  double nu  = MatParams[ 1];
  double m   = MatParams[ 2];
  double c   = MatParams[ 3];
  double phi = MatParams[ 4];
  double Rf  = MatParams[ 5];
  double pa  = MatParams[ 6];
  double rho = MatParams[30];
  double aT  = MatParams[31];
  // Create clone and return
  DuncanChang* newClone =
    new DuncanChang(myID, E, nu, c, phi, m, Rf, pa, rho, aT);
  return newClone;
}
void DuncanChang::set_strain(const Vector& De) {
  eTrial = eTotal+De;
  sTrial=(this->get_C())*eTrial;
}
const Matrix& DuncanChang::get_C() {
  C.clear();
  const Vector& s = sTrial.eigenvalues();
  double s1 = -s[2];
  double s3 = -s[0];
  cout << s1 << '\t' << s3 << endl;

  // Material parameters
  double E   =MatParams[ 0];
  double nu  =MatParams[ 1];
  double m   =MatParams[ 2];
  double c   =MatParams[ 3];
  double phi =MatParams[ 4];
  double Rf  =MatParams[ 5];
  double pa  =MatParams[ 6];
  double d = 1-Rf*(1-sin(phi)*(s1-s3))/(2*c*cos(phi)+2*s3*sin(phi));
  if (s3 < 0) s3 = 0;
  double Et = d*d*E*pa*pow(s3/pa, m);
  if (num::tiny(Et)) Et = E;
  cout << Et << endl;

  // Find and return C
  double Em = Et/((1.+nu)*(1.-2*nu));
  C(0, 0)=Em*(1.-nu);   C(0, 1)=Em*nu;        C(0, 2)=Em*nu;
  C(1, 0)=Em*nu;        C(1, 1)=Em*(1.-nu);   C(1, 2)=Em*nu;
  C(2, 0)=Em*nu;        C(2, 1)=Em*nu;        C(2, 2)=Em*(1.-nu);
  C(3, 3)=Em*0.5*(1.-2*nu);
  C(4, 4)=Em*0.5*(1.-2*nu);
  C(5, 5)=Em*0.5*(1.-2*nu);
  return C;
}
void DuncanChang::commit() {
  eTotal = eTrial;
  sConvg = sTrial;
  this->track();
}
/**
 * Add a record to the tracker.
 * If \a myTracker pointer is null (no tracker is added) just return.
 * Otherwise gather info and send them to the tracker.
 * The domain should be already updated!
 */
void DuncanChang::track() {
  if (myTracker == 0) return;
  // define an output string stream
  std::ostringstream s;
  // start saving
  s << "{";
  // save lambda
  s << "\"lambda\":"  << pD->get_lambda() << ",";
  // save time
  s << "\"time\":"    << pD->get_time_curr() << ",";
  // save self
  s << "\"data\":{";
  s << "\"sigm\":"    << sConvg << ',';
  s << "\"epst\":"    << eTotal << ',';
//  s << "\"epsp\":"    << ePConvg << ',';
  s << "\"epsv\":"    << eTotal[0]+eTotal[1]+eTotal[2] << ',';
  s << "\"p\":"       << sConvg.p() << ',';
  s << "\"q\":"       << sConvg.q();
  s << "}";
  // finalize
  s << "}";
  // convert to c style string and return
  // needs to be converted to a static string before
  /// @todo: check for refactoring
  static string tmp;
  tmp = s.str();
  myTracker->track(tmp.c_str());
}
