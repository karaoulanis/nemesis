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

#include "material/mohr_coulomb.h"
#include <sstream>
#include <string>
#include <vector>
#include "main/nemesis_debug.h"

Matrix MohrCoulomb::C(6, 6, 0.);
Matrix MohrCoulomb::C3(3, 3, 0.);

MohrCoulomb::MohrCoulomb()
    : myElastic(0),
      plastic(false),
      inaccurate(0) {
}

MohrCoulomb::MohrCoulomb(int id, MultiaxialMaterial* elastic, double c,
                         double phi, double alpha, double eta)
    : MultiaxialMaterial(id, 0., 0.),
      plastic(false),
      inaccurate(0) {
  // Get the elastic part
  myElastic = elastic->get_clone();
  MatParams[30] = myElastic->get_param(30);
  MatParams[31] = myElastic->get_param(31);
  // Material properties
  MatParams[0] = c;
  MatParams[1] = phi;
  MatParams[2] = alpha;
  MatParams[3] = eta;
  // Material state
  // ePTrial.Resize(6, 0.); ePConvg.Resize(6, 0.);
  // qTrial.Resize(6, 0.);  qConvg.Resize(6, 0.);
  // aTrial = 0.;           aConvg = 0.;
}


MohrCoulomb::~MohrCoulomb() {
  delete myElastic;
}


MultiaxialMaterial* MohrCoulomb::get_clone() {
  // Material parameters
  int myID     = this->get_id();
  double c     = MatParams[ 0];
  double phi   = MatParams[ 1];
  double alpha = MatParams[ 2];
  double eta   = MatParams[ 3];
  // Create clone and return
  MohrCoulomb* newClone = new MohrCoulomb(myID, myElastic, c, phi, alpha, eta);
  return newClone;
}

/**
 * Update stresses given a total strain increment.
 * @param De Vector containing total strain increment.
 */
void MohrCoulomb::set_strain(const Vector& De, const double Dt) {
  // material properties
  double E  = myElastic->get_param(0);
  double nu = myElastic->get_param(1);
  double c     = MatParams[ 0];
  double phi   = MatParams[ 1];
  double alpha = MatParams[ 2];
  double eta   = MatParams[ 3];

  // derivatives
  std::vector<Vector> df(3);
  df[0].Resize(3);
  df[1].Resize(3);
  df[2].Resize(3);
  // df[3].resize(3);
  df[0][0] = 1.0+sin(phi);
  df[1][0] = 0.0;
  df[2][0] = 1.0+sin(phi);
  // df[3][0]=+1.0;
  df[0][1] = 0.0;
  df[1][1] = 1.0+sin(phi);
  df[2][1] =-1.0+sin(phi);
  // df[3][1]=+1.0;
  df[0][2] =-1.0+sin(phi);
  df[1][2] =-1.0+sin(phi);
  df[2][2] = 0.0;
  // df[3][2]=+1.0;

  std::vector<Vector> dg(3);
  dg[0].Resize(3);
  dg[1].Resize(3);
  dg[2].Resize(3);
  // df[3].resize(3);
  dg[0][0] = 1.0+sin(alpha);
  dg[1][0] = 0.0;
  dg[2][0] = 1.0+sin(alpha);
  // df[3][0]=+1.0;
  dg[0][1] = 0.0;
  dg[1][1] = 1.0+sin(alpha);
  dg[2][1] =-1.0+sin(alpha);
  // df[3][1] =+1.0;
  dg[0][2] =-1.0+sin(alpha);
  dg[1][2] =-1.0+sin(alpha);
  dg[2][2] = 0.0;
  // df[3][2] =+1.0;

  // elasticity matrix
  C3(0, 0) =   1/E;
  C3(0, 1) = -nu/E;
  C3(0, 2) = -nu/E;
  C3(1, 0) = -nu/E;
  C3(1, 1) =   1/E;
  C3(1, 2) = -nu/E;
  C3(2, 0) = -nu/E;
  C3(2, 1) = -nu/E;
  C3(2, 2) =   1/E;

  // spectral decomposition
  static Vector s(3);
  static Matrix sV(3, 3);
  eTrial = eTotal+De;
  sTrial = sConvg+(this->get_C())*De;
  spectralDecomposition(sTrial, s, sV);
  Vector snn = sTrial;

  // yield function
  static Vector f(3);
  f[0] = (s[0]-s[2])+(s[0]+s[2])*sin(phi)-2*c*cos(phi);
  f[1] = (s[1]-s[2])+(s[1]+s[2])*sin(phi)-2*c*cos(phi);
  f[2] = (s[0]-s[1])+(s[0]+s[1])*sin(phi)-2*c*cos(phi);

  // if (f[0]>0. && f[1]>0. && f[2]>0.) f[2]=-1.;

  std::vector<int> active;
  for (unsigned i = 0; i < 3; i++)
  if (f[i]>0.) active.push_back(i);

  // Elastic case
  if (active.empty()) return;
  plastic = true;

  // Plastic case
  static Matrix A;
  static Vector x;
  static Vector R;

  active.clear();
  for (unsigned i = 0; i < 3; i++) {
    if (f[i]>0.) active.push_back(i);
  }

  for (int k = 0; k < 4; k++) {
    A.Resize(3+active.size(), 3+active.size(), 0.);
    x.Resize(3+active.size());
    R.Resize(3+active.size());
    R.Clear();
    A.Append(C3, 0, 0);
    for (unsigned i = 0; i < active.size(); i++) {
      A.AppendCol(dg[active[i]],  0, 3+i);
      A.AppendRow(df[active[i]], 3+i,  0);
      R[3+i]=-f[active[i]];
    }
    // solve
    A.Solve(x, R);
    // check
    bool restart = false;
    for (unsigned i = 0; i < active.size(); i++) {
      if (x[3+i] < 0.) {
        active.erase(active.begin()+i, active.begin()+i+1);
        restart = true;
      }
    }
    if (restart) continue;
    // update
    for (int i = 0; i < 3; i++) s[i]+=x[i];
    break;
  }

  // coordinate transformation
  sTrial[0] = s[0]*sV(0, 0)*sV(0, 0)
             +s[1]*sV(1, 0)*sV(1, 0)
             +s[2]*sV(2, 0)*sV(2, 0);
  sTrial[1] = s[0]*sV(0, 1)*sV(0, 1)
             +s[1]*sV(1, 1)*sV(1, 1)
             +s[2]*sV(2, 1)*sV(2, 1);
  sTrial[2] = s[0]*sV(0, 2)*sV(0, 2)
             +s[1]*sV(1, 2)*sV(1, 2)
             +s[2]*sV(2, 2)*sV(2, 2);
  sTrial[3] = s[0]*sV(0, 0)*sV(0, 1)
             +s[1]*sV(1, 0)*sV(1, 1)
             +s[2]*sV(2, 0)*sV(2, 1);
  sTrial[4] = s[0]*sV(0, 1)*sV(0, 2)
             +s[1]*sV(1, 1)*sV(1, 2)
             +s[2]*sV(2, 1)*sV(2, 2);
  sTrial[5] = s[0]*sV(0, 0)*sV(0, 2)
             +s[1]*sV(1, 0)*sV(1, 2)
             +s[2]*sV(2, 0)*sV(2, 2);

  // check
  f[0]=(s[0]-s[2])+(s[0]+s[2])*sin(phi)-2*c*cos(phi);
  f[1]=(s[1]-s[2])+(s[1]+s[2])*sin(phi)-2*c*cos(phi);
  f[2]=(s[0]-s[1])+(s[0]+s[1])*sin(phi)-2*c*cos(phi);

  // if (f[0]>1e-8 || f[1]>1e-8 || f[2]>1e-8)
  //   cout << "inacurate : "<<s[0]<<"\t"<<s[1]<<"\t"<<s[2]<<"\t"<<endl;
  
  // creep
  if (eta > 0.) {
    sTrial=(snn + Dt / eta * sTrial) / (1. + Dt / eta);
  }
}
/**
 * Commit material state.
 */
void MohrCoulomb::Commit() {
  // report(inaccurate);
  inaccurate = 0;
  eTotal = eTrial;  /// @todo
  sConvg = sTrial;
}
/**
 * Get tangent material matrix.
 * @todo fill it
 * @return A reference to the tangent material matrix.
 */
const Matrix& MohrCoulomb::get_C() {
  return myElastic->get_C();
}

bool MohrCoulomb::isPlastic() {
  return plastic;
}

void MohrCoulomb::Save(std::ostream* os) {
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
