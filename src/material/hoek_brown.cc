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

#include "material/hoek_brown.h"
#include <sstream>
#include <string>
#include "main/nemesis_debug.h"

Matrix HoekBrown::C(6, 6, 0.);
Matrix HoekBrown::C3(3, 3, 0.);

/**
 * Default constructor.
 */
HoekBrown::HoekBrown()
      : myElastic(0),
        plastic(false),
        inaccurate(0),
        f(0),
        dfds(0),
        dgds(0),
        aTrial(0.),
        aConvg(0.) {
}
/**
 * Constructor.
 * @param ID The id of this material
 * @param elasticID The id of the elastic material
 * @param si Uniaxial compressive strength of the impact rock
 * @param sp Empirically determined parameter
 * @param mb Empirically determined parameter
 * @param mbb Empirically determined parameter (??)
 * @param alpha Empirically determined parameter
 */
HoekBrown::HoekBrown(int id, MultiaxialMaterial* elastic, double si, double sp,
                     double mb, double mbb, double alpha, double eta)
    : MultiaxialMaterial(id, 0., 0.),
      plastic(false),
      inaccurate(0),
      aTrial(0.),
      aConvg(0.) {
  // Get the elastic part
  myElastic = elastic->get_clone();
  MatParams[30] = myElastic->get_param(30);
  MatParams[31] = myElastic->get_param(31);

  // Material properties
  MatParams[0] = si;
  MatParams[1] = sp;
  MatParams[2] = mb;
  MatParams[3] = mbb;
  MatParams[4] = alpha;
  MatParams[5] = eta;

  // Material state
  // ePTrial.Resize(6, 0.); ePConvg.Resize(6, 0.);
  // qTrial.Resize(6, 0.);  qConvg.Resize(6, 0.);
  // aTrial = 0.;           aConvg = 0.;

  f.Resize(3, 0.);
  dfds.resize(3);
  dfds[0].Resize(3);
  dfds[1].Resize(3);
  dfds[2].Resize(3);
  dgds.resize(3);
  dgds[0].Resize(3);
  dgds[1].Resize(3);
  dgds[2].Resize(3);
  d2gdsds.resize(3);
  d2gdsds[0].Resize(3, 3, 0.);
  d2gdsds[1].Resize(3, 3, 0.);
  d2gdsds[2].Resize(3, 3, 0.);
}

/**
 * Default destructor.
 * Needs to delete cloned elastic material.
 */
HoekBrown::~HoekBrown() {
  delete myElastic;
}

/**
 * Create a clone of this material.
 */
MultiaxialMaterial* HoekBrown::get_clone() {
  // Material parameters
  int id    = this->get_id();
  double sigma_ci = MatParams[ 0];
  double sp       = MatParams[ 1];
  double mb       = MatParams[ 2];
  double mbb      = MatParams[ 3];
  double alpha    = MatParams[ 4];
  double eta      = MatParams[ 5];
  // Create clone and return
  HoekBrown* clone = new HoekBrown(id, myElastic, sigma_ci, sp, mb, mbb, alpha,
                                   eta);
  return clone;
}

/**
 * Find the value of the gield surfaces w.r.t. to stress.
 * The derivative is stored at Vector g.
 * See http://invenio.lib.auth.gr/record/123972/files/karaoulanis.pdf p.121.
 * @param s Stress Vector.
 */
void HoekBrown::find_f(const Vector& s, double /*q*/) {
  double sigma_ci = MatParams[ 0];
  double sp       = MatParams[ 1];
  double mb       = MatParams[ 2];
  double alpha    = MatParams[ 4];
  f[0]=s[0]-s[2]-sigma_ci*pow(sp-mb*s[0]/sigma_ci, alpha);
  f[1]=s[1]-s[2]-sigma_ci*pow(sp-mb*s[1]/sigma_ci, alpha);
  f[2]=s[0]-s[1]-sigma_ci*pow(sp-mb*s[0]/sigma_ci, alpha);
}

/**
 * Find the first derivative of the plastic potential surfaces w.r.t. to stress.
 * The derivative is stored at std::vector<Vector> dfds.
 * See http://invenio.lib.auth.gr/record/123972/files/karaoulanis.pdf p.121.
 * @param s Stress Vector.
 */
void HoekBrown::find_dfds(const Vector& s, double /*q*/) {
  double sigma_ci = MatParams[ 0];
  double sp       = MatParams[ 1];
  double mb       = MatParams[ 2];
  double alpha    = MatParams[ 4];
  double c0 = alpha*mb*pow(sp-mb*s[0]/sigma_ci, alpha-1);
  double c1 = alpha*mb*pow(sp-mb*s[1]/sigma_ci, alpha-1);
  dfds[0][0] =  1.0+c0;
  dfds[1][0] =  0.0;
  dfds[2][0] =  1.0+c0;
  dfds[0][1] =  0.0;
  dfds[1][1] =  1.0+c1;
  dfds[2][1] = -1.0;
  dfds[0][2] = -1.0;
  dfds[1][2] = -1.0;
  dfds[2][2] =  0.0;
}

/**
 * Find the first derivative of the plastic potential surfaces w.r.t. to stress.
 * The derivative is stored at std::vector<Vector> dgds.
 * See http://invenio.lib.auth.gr/record/123972/files/karaoulanis.pdf p.121.
 * @param s Stress Vector.
 */
void HoekBrown::find_dgds(const Vector& s, double /*q*/) {
  double sigma_ci = MatParams[ 0];
  double sp       = MatParams[ 1];
  double mb       = MatParams[ 2];
  double alpha    = MatParams[ 4];
  double c0 = alpha*mb*pow(sp-mb*s[0]/sigma_ci, alpha-1);
  double c1 = alpha*mb*pow(sp-mb*s[1]/sigma_ci, alpha-1);
  dgds[0][0]= 1.0+c0;
  dgds[1][0]= 0.0;
  dgds[2][0]= 1.0+c0;
  dgds[0][1]= 0.0;
  dgds[1][1]= 1.0+c1;
  dgds[2][1]=-1.0;
  dgds[0][2]=-1.0;
  dgds[1][2]=-1.0;
  dgds[2][2]= 0.0;
}

/**
 * Find the second derivative of the yield  potential surfaces w.r.t. to stress.
 * The derivative is stored at std::vector<Matrix> d2gdsds.
 * See http://invenio.lib.auth.gr/record/123972/files/karaoulanis.pdf p.121-122.
 * @param s Stress Vector.
 */
void HoekBrown::find_d2gdsds(const Vector& s, double /*q*/) {
  double sigma_ci = MatParams[ 0];
  double sp       = MatParams[ 1];
  double mb       = MatParams[ 2];
  double alpha    = MatParams[ 4];
  double c00 = alpha*(1-alpha)*mb*mb*pow(sp-mb*s[0]/sigma_ci, alpha-2)/sigma_ci;
  double c11 = alpha*(1-alpha)*mb*mb*pow(sp-mb*s[1]/sigma_ci, alpha-2)/sigma_ci;
  d2gdsds[0](0, 0)=c00;
  d2gdsds[1](1, 1)=c11;
  d2gdsds[2](0, 0)=c00;
}

/**
 * Update stresses given a total strain increment.
 * @param De Vector containing total strain increment.
 */
void HoekBrown::set_strain(const Vector& De, const double Dt) {
  // material properties
  double E = myElastic->get_param(0);
  double nu= myElastic->get_param(1);

  // elasticity matrix
  C3(0, 0)=  1/E;
  C3(0, 1)=-nu/E;
  C3(0, 2)=-nu/E;
  C3(1, 0)=-nu/E;
  C3(1, 1)=  1/E;
  C3(1, 2)=-nu/E;
  C3(2, 0)=-nu/E;
  C3(2, 1)=-nu/E;
  C3(2, 2)=  1/E;

  // spectral decomposition
  Vector s(3), De3(3);
  Vector sTrial3(3);  /// @todo: remove
  Matrix sV(3, 3), eV(3, 3);
  aTrial = aConvg;
  eTrial = eTotal+De;
  sTrial = sConvg+(this->get_C())*De;
  spectralDecomposition(sTrial, s, sV);
  Vector eTrial3 = C3*s;
  sTrial3 = s;
  Vector snn = sTrial;

  // Yield functions
  double q = 0.;
  this->find_f(s, q);
  // report(f, "f", 8, 4);

  // Setup
  int nf = f.get_size();
  std::vector<bool> active(nf);
  std::vector<double> DLambda(nf, 0.);
  Matrix A;
  Vector R;
  Vector x;

  // Find active surfaces
  int nActive = 0;
  for (int i = 0; i < nf; i++) {
    if (f[i]>0.) {
      active[i]=true;
      nActive++;
    } else {
      active[i]=false;
    }
  }

  // Check for purely elastic
  if (nActive == 0) {
    return;
  }

  plastic = true;

  // Check for stupid
  if (nActive == 3) {
    // return;
  }

  // Now its time for plasticity
  this->find_dfds(s, q);
  this->find_dgds(s, q);
  this->find_d2gdsds(s, q);

  for (; ; ) {
    // setup jacobian
    A.Resize(3+nActive, 3+nActive, 0.);
    R.Resize(3+nActive, 0.);
    x.Resize(3+nActive);
    A.Append(C3, 0, 0);
    R.Append(-C3*s+eTrial3, 0, 1., 0.);
    int pos = 0;
    for (int i = 0; i < 3; i++) {
      if (active[i]) {
        A.Append(DLambda[i]*d2gdsds[i],    0,     0, 1., 1.);
        A.AppendCol(dgds[i],               0, 3+pos, 1., 0.);
        A.AppendRow(dfds[i],           3+pos,     0, 1., 0.);
        R.Append(-DLambda[i]*dgds[i],             0, 1., 1.);
        R[3+pos]=-f[i];
        pos++;
      }
    }

    // check for convergence
    if (R.Twonorm() < 1.e-09) {
      bool restart = false;
      pos = 0;
      for (int i = 0; i < 3; i++) {
        if (active[i]) {
          if (DLambda[i] < 0.) {
            active[i]=false;
            nActive--;
            restart = true;
          }
          pos++;
        }
      }
      if (restart) {
        // initialize stresses and lambdas
        s = sTrial3;
        for (int i = 0; i < nf; i++) DLambda[i]=0.;
        this->find_f(s, q);
        this->find_dfds(s, q);
        this->find_dgds(s, q);
        restart = false;
        // restart with initial data
        continue;
      } else {
        break;
      }
    }

    // solve
    A.Solve(x, R);

    // Update stresses
    for (int i = 0; i < 3; i++)
      s[i]+=x[i];
    // Update dLambda
    pos = 0;
    for (int i = 0; i < 3; i++) {
      if (active[i]) {
        DLambda[i]+=x[3+pos];
        pos++;
      }
    }
    // update vectors
    this->find_f(s, q);
    this->find_dfds(s, q);
    this->find_dgds(s, q);
    this->find_d2gdsds(s, q);
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
  // creep
  // double Dt  = 10000.;
  // double eta = 100000.;
    //double Dt  = 1.;
 
    
    //sTrial=(snn+(Dt/eta)*sTrial)/(1+Dt/eta);
  // std::cout << sTrial[1] << std::endl;
  //aTrial=(ann+(Dt/eta)*aTrial)/(1+Dt/eta);
  //double Dt_eta = 0.02;
  double eta = MatParams[5];
  if (eta > 0.) {
    sTrial=(snn + Dt / eta * sTrial) / (1. + Dt / eta);
  }
}

/**
 * Commit material state.
 */
void HoekBrown::Commit() {
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
const Matrix& HoekBrown::get_C() {
  return myElastic->get_C();
}

bool HoekBrown::isPlastic() {
  return plastic;
}

void HoekBrown::Save(std::ostream* os) {
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
