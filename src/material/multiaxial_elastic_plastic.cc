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

#include "material/multiaxial_elastic_plastic.h"
#include <sstream>
#include <string>
#include "containers/containers.h"
#include "main/nemesis_debug.h"
#include "material/evolution_law.h"
#include "material/linear_equivalent_el.h"
#include "material/surface.h"

Matrix MultiaxialElastoPlastic::C(6, 6, 0.);

MultiaxialElastoPlastic::MultiaxialElastoPlastic()
    : myElastic(0),
      ePTrial(0),
      ePConvg(0),
      qTrial(0),
      qConvg(0),
      aTrial(0.),
      aConvg(0.),
      plastic(false),
      nHardeningVariables(0),
      fSurfaces(0),
      gSurfaces(0),
      EL(0) {
}

MultiaxialElastoPlastic::
MultiaxialElastoPlastic(int id, MultiaxialMaterial* elastic)
    : MultiaxialMaterial(id, 0., 0.),
      aTrial(0.),
      aConvg(0.),
      plastic(false),
      nHardeningVariables(0),
      fSurfaces(0),
      gSurfaces(0),
      EL(new LinearEquivalentEL()) {
  // Get the elastic part
  myElastic = elastic->get_clone();
  MatParams[30] = myElastic->get_param(30);
  MatParams[31] = myElastic->get_param(31);
  // Material state
  ePTrial.Resize(6, 0.);
  ePConvg.Resize(6, 0.);
  qTrial.Resize(6, 0.);
  qConvg.Resize(6, 0.);
}

MultiaxialElastoPlastic::~MultiaxialElastoPlastic() {
  delete myElastic;
  Containers::vector_delete(&fSurfaces);
  Containers::vector_delete(&gSurfaces);
}
/**
 * Update stresses given a total strain increment.
 * @param De Vector containing total strain increment.
 */
void MultiaxialElastoPlastic::set_strain(const Vector& De) {
  // this->returnMapTest(De);
  this->returnMapSYS2(De);
  // this->returnMapMYS3(De);
  // if (fSurfaces.size()==1) this->returnMapSYS(De);
  // if (fSurfaces.size()==1) this->returnMapTest(De);
  // else          this->returnMapMYS(De);
}
/**
 * Commit material state.
 */
void MultiaxialElastoPlastic::Commit() {
  eTotal = eTrial;  /// @todo
  ePConvg = ePTrial;
  sConvg = sTrial;
  qConvg = qTrial;
  aConvg = aTrial;
}

void MultiaxialElastoPlastic::returnMapTest(const Vector& De) {
  int nIter = 20;
  double tol1 = 1e-6;
  double tol2 = 1e-6;
  Vector R(7, 0.);
  Matrix Cel = myElastic->get_C();
  Matrix invCel = Inverse(Cel);
  Matrix A(7, 7, 0.);
  ePTrial = ePConvg;
  aTrial = aConvg;
  double dg = 0;
  double q = 500.*aTrial;
  Surface* fS = fSurfaces[0];
  Surface* gS = gSurfaces[0];
  // double eta = 1000.;
  // double dt = pD->get_time_incr();

  // ===========================================================================
  // Step 1: Compute trial stress
  // ===========================================================================
  eTrial = eTotal+De;
  sTrial = Cel*(eTrial-ePTrial);

  // ===========================================================================
  // Step 2: Check for plastic process
  // ===========================================================================
  if (fS->get_f(sTrial, q) <= 0) return;
  for (int k = 0; k < nIter; k++) {
    q = 500.*aTrial;
    // Residual Vector
    Vector tmp(6, 0.);
    tmp = De;
    tmp-=invCel*(sTrial-sConvg);
    tmp-=dg*(gS->get_dfds(sTrial, aTrial));
    R.Clear();
    R.Append(tmp, 0);
    R[6]=-(fS->get_f(sTrial, q));
    // Convergence check
    if ((tmp.Twonorm() < tol1) && (std::abs(fS->get_f(sTrial, q)) < tol2)) {
      aTrial+=dg;
      ePTrial+=dg*gS->get_dfds(sTrial, aTrial);
      std::cout << dg << std::endl;
      break;
    }
    // Jacobian
    A.Clear();
    A.Append(invCel+dg*(gS->get_df2dss(sTrial, aTrial)), 0, 0);
    for (int i = 0;i < 6;i++) A(6, i)=(fS->get_dfds(sTrial, aTrial))[i];
    for (int i = 0;i < 6;i++) A(i, 6)=(gS->get_dfds(sTrial, aTrial))[i];
    // efac
    Vector r = fS->get_dfds(sTrial, aTrial);
    double rm = r[0]+r[1]+r[2];
    r[0]-=rm;
    r[1]-=rm;
    r[2]-=rm;
    double efac = sqrt(2./3.)*sqrt(r[0]*r[0]
                                  +r[1]*r[1]
                                  +r[2]*r[2]
                                  +0.5*(r[3]*r[3]+r[4]*r[4]+r[5]*r[5]));
    // H
    double H=-500.;
    A(6, 6)=H*efac;
    Vector x(7, 0.);
    A.Solve(x, R);
    for (int i = 0;i < 6;i++) sTrial[i]+=x[i];
    dg+=x[6];
    aTrial+=dg*efac;
  }
  // sTrial=(ss+(dt/eta)*sTrial)/(1+(dt/eta));
}

void MultiaxialElastoPlastic::returnMapSYS(const Vector& De) {
  // ===========================================================================
  // Setup
  // ===========================================================================
  int nIter = 25;
  double tol1 = 1e-6;
  double tol2 = 1e-6;
  Matrix Cel  = myElastic->get_C();
  Matrix invCel = Inverse(Cel);

  aTrial  = aConvg;
  ePTrial = ePConvg;
  eTrial  = eTotal+De;
  double dg = 0.;

  Surface* fS = fSurfaces[0];
  Surface* gS = gSurfaces[0];

  // ===========================================================================
  // Step 1: Compute trial stress
  // ===========================================================================
  sTrial = Cel*(eTrial-ePTrial);

  // ===========================================================================
  // Step 2: Check for plastic process
  // ===========================================================================
  if (fS->get_f(sTrial, aTrial) <=0 ) return;

  for (int k = 0; k < nIter; k++) {
    // =========================================================================
    // Step 3: Evaluate flow rule, hardening law and residuals
    // =========================================================================
    sTrial = Cel*(eTrial-ePTrial);

    static Vector R(6, 0.);
    if (nHardeningVariables>0)
      R.Resize(7, 0.);
    R.Clear();
    Vector temp = De;
    temp-=invCel*(sTrial-sConvg);
    temp-=dg*(gS->get_dfds(sTrial, aTrial));
    R.Append(temp, 0);
    if (nHardeningVariables>0) {
      R[6]=-aTrial+aConvg+dg*(EL->get_h(gS->get_dfds(sTrial, aTrial)));
    }
    // =========================================================================
    // Step 4: Check convergence
    // =========================================================================
    if (R.Twonorm() < tol1 && std::abs(fS->get_f(sTrial, aTrial)) < tol2) break;

    // =========================================================================
    // Step 5: Compute elastic moduli and consistent tangent moduli
    // =========================================================================
    // Matrix A
    Matrix A(6+nHardeningVariables, 6+nHardeningVariables, 0.);
    A.Append(invCel+dg*(gS->get_df2dss(sTrial, aTrial)), 0, 0);
    if (nHardeningVariables>0) {
      Vector v = dg*(gS->get_df2dsa(sTrial, aTrial));
      for (int i = 0;i < 6;i++) A(i, 6)=v[i];
      for (int i = 0;i < 6;i++) A(6, i)=v[i];
      A(6, 6)=-1.+EL->get_dhda(sTrial, ePTrial);
    }
    // A.report();
    A = Inverse(A);

    // =========================================================================
    // Step 6: Obtain increment to consistency parameter
    // =========================================================================
    Vector vf(6+nHardeningVariables, 0);
    vf.Append(fS->get_dfds(sTrial, aTrial), 0);
    for (int i = 0;i < nHardeningVariables;i++)
      vf[6+i]=fS->get_dfda(sTrial, aTrial);

    Vector vh(6+nHardeningVariables, 0);
    vh.Append(fS->get_dfds(sTrial, aTrial), 0);
    for (int i = 0;i < nHardeningVariables;i++) {
      vh[6+i]=EL->get_h(gS->get_dfds(sTrial, aTrial));
    }
    double fy = fS->get_f(sTrial, aTrial);
    double ddg=(fy-vf*(A*R))/(vf*(A*vh));

    // =========================================================================
    // Step 7: Obtain incremental plastic strains and internal variables
    // =========================================================================
    Vector x(6+nHardeningVariables, 0.);
    x = A*(-R-ddg*vh);

    Vector dEp(6, 0.);
    for (int i = 0;i < 6;i++) dEp[i]=x[i];
    dEp=-invCel*dEp;
        double da = x[6];

    // =========================================================================
    // Step 8: Update
    // =========================================================================
    ePTrial+=dEp;
    aTrial+=da;
    dg+=ddg;
  }
  // std::cout << aTrial << '\t' << dg << '\t'
  //      << (r[0]*r[0]+r[1]*r[1]+r[2]*r[2]
  //          +0.5*r[3]*r[3]+0.5*r[4]*r[4]+0.5*r[5]*r[5]) << std::endl;
}
double kappa(const Vector& v, double Dlambda) {
  double eq = 2./3.*(v[0]*v[0]+    v[1]*v[1]+    v[2]*v[2]
          +0.5*v[3]*v[3]+0.5*v[4]*v[4]+0.5*v[5]*v[5]);
  // std::cout << eq << std::endl;
  return Dlambda*sqrt(eq);
}
double dkds(const Vector& /*s*/, double /*Dlambda*/) {
  return 0.;
}
double dkdl(const Vector& v, double /*Dlambda*/) {
  double eq = 2./3.*(v[0]*v[0]+    v[1]*v[1]+    v[2]*v[2]
          +0.5*v[3]*v[3]+0.5*v[4]*v[4]+0.5*v[5]*v[5]);
  // std::cout << eq << std::endl;
  return sqrt(eq);
}
void MultiaxialElastoPlastic::returnMapSYS2(const Vector& De) {
  // ===========================================================================
  // Setup
  // ===========================================================================
  int nIter = 25;
  double tol1 = 1e-6;
  double tol2 = 1e-6;
  Matrix Cel    = myElastic->get_C();
  Matrix invCel = Inverse(Cel);

  aTrial  = aConvg;
  ePTrial = ePConvg;
  eTrial  = eTotal+De;
  double dg = 0.;

  Surface* fS = fSurfaces[0];
  Surface* gS = gSurfaces[0];

  //=========================================================================
  // Step 1: Compute trial stress
  //=========================================================================
  sTrial = Cel*(eTrial-ePTrial);
  Vector ss(6, 0.);
  ss = sTrial;
  //=========================================================================
  // Step 2: Check for plastic process
  //=========================================================================
  if (fS->get_f(sTrial, aTrial) <= 0) return;

  for (int k = 0; k < nIter; k++) {
    //=====================================================================
    // Step 3: Evaluate flow rule, hardening law and residuals
    //=====================================================================
    sTrial = Cel*(eTrial-ePTrial);

    static Vector R(8, 0.);
    Vector temp = De;
    temp-=invCel*(sTrial-sConvg);
    temp-=dg*(gS->get_dfds(sTrial, aTrial));
    R.Append(temp, 0);
    R[6]=-aTrial+aConvg+kappa(gS->get_dfds(sTrial, aTrial), dg);
    R[7]=-(fS->get_f(sTrial, aTrial));

    //=====================================================================
    // Step 4: Check convergence
    //=====================================================================
    if (R.Twonorm() < tol1 && std::abs(fS->get_f(sTrial, aTrial)) < tol2) break;

    //=====================================================================
    // Step 5: Compute elastic moduli and consistent tangent moduli
    //=====================================================================
    // Matrix A
    Matrix A(8, 8, 0.);
    A.Append(invCel+dg*(gS->get_df2dss(sTrial, aTrial)), 0, 0);
    for (int i = 0;i < 6;i++) A(i, 6)=dg*(gS->get_df2dsa(sTrial, aTrial)[i]);
    for (int i = 0;i < 6;i++) A(i, 7)=gS->get_dfds(sTrial, aTrial)[i];
    for (int i = 0;i < 6;i++) A(6, i)=0.;
    A(6, 6)=-1.;
    A(6, 7)=dkdl(gS->get_dfds(sTrial, aTrial), dg);
    for (int i = 0;i < 6;i++) A(7, i)=fS->get_dfds(sTrial, aTrial)[i];
    A(7, 6)=fS->get_dfda(sTrial, aTrial);
    A(7, 7)=0.;

    // A.report();
    Vector x(8, 0.);
    A.Solve(x, R);

    Vector dEp(6, 0.);
    for (int i = 0;i < 6;i++) dEp[i]=x[i];
    dEp=-invCel*dEp;
        double da = x[6];
    double ddg = x[7];

    //=====================================================================
    // Step 8: Update
    //=====================================================================
    ePTrial+=dEp;
    aTrial+=da;
    dg+=ddg;
  }
//  double eta = 1000.;
//  double Dt = pD->get_time_incr();
//  sTrial=(ss    +(Dt/eta)*sTrial)/(1+Dt/eta);
//  aTrial=(aConvg+(Dt/eta)*aTrial)/(1+Dt/eta);
//  ePTrial = eTrial-invCel*sTrial;
}

/**
 * Multiple surface return mapping.
 * @param De Vector containing total strain increment.
 */
void MultiaxialElastoPlastic::returnMapMYS(const Vector& De) {
  //=========================================================================
  // Setup
  //=========================================================================
  int nIter = 15;
  double tol1 = 1e-6;
  double tol2 = 1e-6;
  Vector R(6, 0.);
  Vector dEp(6, 0.);
  Matrix Cel   =myElastic->get_C();
  Matrix invCel = Inverse(Cel);
  Vector ddg(12, 0.);
  Vector dg(12, 0.);
  double q = 0.;

  ePTrial = ePConvg;   // todo:CHECK!!!!
  // aTrial = aConvg;  // todo:CHECK!!!!
  // aTrial.Clear();   // todo:CHECK!!!!
  eTotal+=De;
  Vector enn = invCel*sConvg+ePConvg+De;
  // Vector enn = eTotal;
  // ===========================================================================
  // Step 1: Compute trial stress
  // ===========================================================================
  // sTrial = sConvg+Cel*De;
  sTrial = Cel*(enn-ePTrial);
  Vector ss = sTrial;
  plastic = false;

  // ===========================================================================
  // Step 2: Check for plastic process
  // ===========================================================================
  int nActiveSurfaces = 0;
  for (unsigned i = 0;i < fSurfaces.size();i++) {
    if (fSurfaces[i]->get_f(sTrial, q)>0) {
      fSurfaces[i]->set_active(true);
      nActiveSurfaces++;
    }
  }
  if (nActiveSurfaces == 0) return;
  // std::cout << nActiveSurfaces << std::endl;

  plastic = true;
  std::cout << std::endl;
  for (int k = 0; k < nIter; k++) {
    // =========================================================================
    // Step 3: Evaluate flow rule, hardening law and residuals
    // =========================================================================
    sTrial = Cel*(enn-ePTrial);
    R=-ePTrial+ePConvg;
    for (unsigned b = 0;b < fSurfaces.size();b++)
      if (fSurfaces[b]->isActive())
        R+=dg[b]*(fSurfaces[b]->get_dfds(sTrial, aTrial));

    nActiveSurfaces = 0;
    for (unsigned b = 0;b < fSurfaces.size();b++) {
      if (fSurfaces[b]->isActive()) {
        nActiveSurfaces++;
      }
    }
    std::cout << "active : " << nActiveSurfaces << std::endl;

    // =========================================================================
    // Step 4: Check convergence
    // =========================================================================
    bool converged = false;
    if (R.Twonorm()>tol1) converged = false;
    for (unsigned a = 0;a < fSurfaces.size();a++) {
      if ((fSurfaces[a]->isActive()) &&
          (std::abs(fSurfaces[a]->get_f(sTrial, q)) < tol2)) {
        converged = true;
      }
    }
    if (converged) break;

    // =========================================================================
    // Step 5: Compute elastic moduli and consistent tangent moduli
    // =========================================================================
    // Matrix A
    Matrix A(6, 6, 0.);
    A.Append(invCel, 0, 0);
    for (unsigned b = 0;b < fSurfaces.size();b++)
      if (fSurfaces[b]->isActive())
        A+=dg[b]*(fSurfaces[b]->get_df2dss(sTrial, aTrial));
    A = Inverse(A);

    // Matrix G
    Matrix Gab(4, 4, 0.);
    Gab(0, 0)=1.0;
    Gab(1, 1)=1.0;
    Gab(2, 2)=1.0;
    Gab(3, 3)=1.0;
    for (unsigned a = 0;a < fSurfaces.size();a++) {
      for (unsigned b = 0;b < fSurfaces.size();b++) {
        if (fSurfaces[a]->isActive() && fSurfaces[b]->isActive()) {
          Vector va(6, 0.);
          Vector vb(6, 0.);
          va = fSurfaces[a]->get_dfds(sTrial, aTrial);
          vb = fSurfaces[b]->get_dfds(sTrial, aTrial);
                  Gab(a, b)=(va*(A*vb));
          // cout << sTrial << std::endl;
          // cout << va << std::endl;
          // cout << vb << std::endl;
          // cout << A << std::endl;
        }
      }
    }
    // std::cout << Gab << std::endl;
    Gab = Inverse(Gab);

    // =========================================================================
    // Step 6: Obtain increment to consistency parameter
    // =========================================================================
    ddg.Clear();
    for (unsigned a = 0;a < fSurfaces.size();a++) {
      for (unsigned b = 0;b < fSurfaces.size();b++) {
        if (fSurfaces[a]->isActive() && fSurfaces[b]->isActive()) {
          Vector vb = fSurfaces[b]->get_dfds(sTrial, aTrial);
          double fb = fSurfaces[b]->get_f(sTrial, q);
          ddg[a]+=Gab(a, b)*(fb-vb*(A*R));
        }
      }
    }

    bool reset = false;
    for (unsigned a = 0; a < fSurfaces.size(); a++) {
      // cout << "a = "<<a<<"  "<<dg[a]+ddg[a]<<std::endl;
      if (fSurfaces[a]->isActive() &&  (dg[a]+ddg[a] < 0.)) {
        fSurfaces[a]->set_active(false);
        dg[a]=0.;
        ddg[a]=0.;
        reset = true;
      }
    }
    if (reset) continue;
    std::cout << k << " dg " << ddg << std::endl;

    // =========================================================================
    // Step 7: Obtain incremental plastic strains and internal variables
    // =========================================================================
    static Vector temp(6);
    temp.Clear();
    for (unsigned b = 0;b < fSurfaces.size();b++)
      if (fSurfaces[b]->isActive())
        temp+=ddg[b]*(fSurfaces[b]->get_dfds(sTrial, aTrial));

    // =========================================================================
    // Step 8: Update
    // =========================================================================
    dEp = invCel*A*(R+temp);
    ePTrial+=dEp;
    for (unsigned a = 0;a < fSurfaces.size();a++)
      if (fSurfaces[a]->isActive()) dg[a]+=ddg[a];
  }
  // if (k == nIter) cout << "FAILED"<<std::endl;
  // if (k == nIter) {
  //   cout.precision(12);
  //   cout << sConvg+Cel*De << std::endl;
  //   cout << De << std::endl;
  //   cout << sConvg << std::endl;
  //   cout << ePConvg << std::endl;
  //   cout << ePTrial << std::endl;
  //   cout << "PASSED : "<<nActiveSurfaces << std::endl;
  // }
  // cout << std::endl;:
  // double dt = pD->get_time_incr();
  // double eta = 1000.;
  // sTrial=((sConvg+Cel*De)+(dt/eta)*(sConvg+Cel*De))/(1+(dt/eta));
  // sTrial=(ss+(dt/eta)*sTrial)/(1+(dt/eta));
}
/**
 * Get tangent material matrix.
 * @todo fill it
 * @return A reference to the tangent material matrix.
 */
const Matrix& MultiaxialElastoPlastic::get_C() {
  return myElastic->get_C();
}

void MultiaxialElastoPlastic::Save(std::ostream* os) {
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

void MultiaxialElastoPlastic::returnMapMYS2(const Vector& De) {
  // ===========================================================================
  // Setup
  // ===========================================================================
  int nIter = 25;
  double tol1 = 1e-6;
  double tol2 = 1e-6;
  Vector R(6, 0.);
  Vector dEp(6, 0.);
  Matrix Cel   =myElastic->get_C();
  Matrix invCel = Inverse(Cel);
  Vector dg(12, 0.);
  double q = 0;

  ePTrial = ePConvg;   // todo:CHECK!!!!
  // aTrial = aConvg;  // todo:CHECK!!!!
  // aTrial.Clear();   // todo:CHECK!!!!
  eTrial = eTotal;
  eTrial+=De;
  Vector enn = invCel*sConvg+ePConvg+De;
  // ===========================================================================
  // Step 1: Compute trial stress
  // ===========================================================================
  // sTrial = sConvg+Cel*De;
  sTrial = Cel*(enn-ePTrial);
  Vector ss = sTrial;
  plastic = false;

  // ===========================================================================
  // Step 2: Check for plastic process
  // ===========================================================================
  int nActiveSurfaces = 0;
  for (unsigned i = 0; i < fSurfaces.size(); i++) {
    if (fSurfaces[i]->get_f(sTrial, q)>0) {
      fSurfaces[i]->set_active(true);
      nActiveSurfaces++;
    }
  }
  if (nActiveSurfaces == 0) return;
  // cout << nActiveSurfaces << std::endl;
  // cout << enn-eTrial << std::endl; /// @todo Check out which is better.

  plastic = true;
  // cout << std::endl;
  for (int k = 0; k < nIter; k++) {
    //=====================================================================
    // Step 3: Evaluate flow rule, hardening law and residuals
    //=====================================================================
    sTrial = Cel*(enn-ePTrial);
    R=-ePTrial+ePConvg;
    for (unsigned b = 0;b < fSurfaces.size();b++)
      if (fSurfaces[b]->isActive())
        R+=dg[b]*(fSurfaces[b]->get_dfds(sTrial, aTrial));

    nActiveSurfaces = 0;
    for (unsigned b = 0;b < fSurfaces.size();b++) {
      if (fSurfaces[b]->isActive()) {
        nActiveSurfaces++;
      }
    }
    // cout << "active : " << nActiveSurfaces << std::endl;

    //==========================================================================
    // Step 4: Check convergence
    //==========================================================================
    bool converged = false;
    if (R.Twonorm()>tol1) converged = false;
    for (unsigned a = 0;a < fSurfaces.size();a++) {
      if ((fSurfaces[a]->isActive()) &&
          (std::abs(fSurfaces[a]->get_f(sTrial, q)) < tol2)) {
        converged = true;
      }
    }
    if (converged) break;
    // cout << "Res:"<<R.Twonorm()<<' '<<abs(fSurfaces[0]->get_f(sTrial))<<endl;

    // =========================================================================
    // Step 5: Compute elastic moduli and consistent tangent moduli
    // =========================================================================
    // Matrix A
    Matrix A(6, 6, 0.);
    A.Append(invCel, 0, 0);
    for (unsigned b = 0;b < fSurfaces.size();b++)
      if (fSurfaces[b]->isActive())
        A+=dg[b]*(fSurfaces[b]->get_df2dss(sTrial, aTrial));
    A = Inverse(A);

    // Matrix G
    Matrix Gab(nActiveSurfaces, nActiveSurfaces, 0.);
    for (unsigned a = 0;a < fSurfaces.size();a++)
      for (unsigned b = 0;b < fSurfaces.size();b++)
        if (fSurfaces[a]->isActive() && fSurfaces[b]->isActive()) {
          Vector va(6, 0.);
          Vector vb(6, 0.);
          va = fSurfaces[a]->get_dfds(sTrial, aTrial);
          vb = fSurfaces[b]->get_dfds(sTrial, aTrial);
                  Gab(a, b)=(va*(A*vb));
          // cout << a << ' ' << b << std::endl;
          // cout << sTrial << std::endl;
          // cout << va << std::endl;
          // cout << vb << std::endl;
          // cout << A << std::endl;
        }

    //=====================================================================
    // Step 6: Obtain increment to consistency parameter
    //=====================================================================
    Vector ddg(nActiveSurfaces, 0.);
    Vector rhs(nActiveSurfaces, 0.);
    for (unsigned a = 0;a < fSurfaces.size();a++) {
      for (unsigned b = 0;b < fSurfaces.size();b++) {
        if (fSurfaces[a]->isActive() && fSurfaces[b]->isActive()) {
          Vector vb = fSurfaces[b]->get_dfds(sTrial, aTrial);
          double fb = fSurfaces[b]->get_f(sTrial, q);
          rhs[a]+=(fb-vb*(A*R));
        }
      }
    }
    bool reset = false;
    if (std::abs(Det(Gab)) < 1e-8) {
      fSurfaces[nActiveSurfaces-1]->set_active(false);
      reset = true;
    }
    // cout << k << " Gab "<< Gab << std::endl;
    // cout << k << " rhs "<< rhs << std::endl;
    // cout << "det " << det(Gab) << std::endl;
    Gab.Solve(ddg, rhs);
    // cout << k<< "ddg " << ddg << std::endl;

    for (unsigned a = 0; a < fSurfaces.size(); a++) {
      // cout << "a = " << a << "  " << dg[a]+ddg[a] << std::endl;
      if (fSurfaces[a]->isActive() && (dg[a]+ddg[a] < 0.)) {
        fSurfaces[a]->set_active(false);
        dg[a]  = 0.;
        ddg[a] = 0.;
        reset  = true;
      }
    }
    if (reset) continue;

    // =========================================================================
    // Step 7: Obtain incremental plastic strains and internal variables
    // =========================================================================
    static Vector temp(6);
    temp.Clear();
    for (unsigned b = 0;b < fSurfaces.size();b++)
      if (fSurfaces[b]->isActive())
        temp+=ddg[b]*(fSurfaces[b]->get_dfds(sTrial, aTrial));

    // =========================================================================
    // Step 8: Update
    // =========================================================================
    dEp = invCel*A*(R+temp);
    ePTrial+=dEp;
    for (unsigned a = 0;a < fSurfaces.size();a++)
      if (fSurfaces[a]->isActive()) dg[a]+=ddg[a];
  }
  // if (k == nIter) cout << "FAILED"<<std::endl;
}
void MultiaxialElastoPlastic::updateStateVariable() {
  // double eq = sqrt(2/3.)*ePTrial.Twonorm();
}

void MultiaxialElastoPlastic::returnMapMYS3(const Vector& De) {
  //=========================================================================
  // Setup
  //=========================================================================
  int nIter = 25;
  double tol1 = 1e-6;
  double tol2 = 1e-9;
  Vector dEp(6, 0.);
  Matrix Cel   =myElastic->get_C();
  Matrix invCel = Inverse(Cel);

  aTrial =aConvg;
  ePTrial = ePConvg;
  eTrial = eTotal+De;
  Vector dg(12, 0.);

  //=========================================================================
  // Step 1: Compute trial stress
  //=========================================================================
  sTrial = sConvg+Cel*De;
  Vector ss(6, 0.);
  ss = sTrial;

  //=========================================================================
  // Step 2: Check for plastic process
  //=========================================================================
  std::vector<int> activeSurfaces;
  activeSurfaces.resize(0);
  for (unsigned i = 0;i < fSurfaces.size();i++)
    if (fSurfaces[i]->get_f(sTrial, aTrial)>tol2) {
      fSurfaces[i]->set_active(true);
      activeSurfaces.push_back(i);
    }
  if (activeSurfaces.empty()) return;
  plastic = true;
  int k = 0;
  for (k = 0; k < nIter; k++) {
    // cout << k<<' '<<activeSurfaces.size()<<std::endl;
    // =========================================================================
    //  Step 3: Evaluate flow rule, hardening law and residuals
    // =========================================================================
    // sTrial = Cel*(eTrial-ePTrial);
    activeSurfaces.resize(0);
    for (unsigned i = 0; i < fSurfaces.size(); i++) {
      if (fSurfaces[i]->get_f(sTrial, aTrial)>tol2) {
        fSurfaces[i]->set_active(true);
        activeSurfaces.push_back(i);
      }
    }

    if (activeSurfaces.size()>tol2 &&
        std::abs(sTrial.Theta())*180./num::pi>29.99) {
      activeSurfaces.resize(0);
      for (unsigned i = 0;i < fSurfaces.size(); i++)
        fSurfaces[i]->set_active(false);
      fSurfaces[0]->set_active(true);
      activeSurfaces.push_back(0);
    }

    // cout << sTrial.theta()/num::pi*180. << std::endl;
    // cout << "active surfaces : "<<activeSurfaces.size()<<std::endl;

    Vector R(7+activeSurfaces.size(), 0.);
    // ----- R[0-5]
    Vector temp = De;
    temp-=invCel*(sTrial-sConvg);
    for (unsigned b = 0;b < fSurfaces.size();b++)
      if (fSurfaces[b]->isActive())
        temp-=dg[b]*(fSurfaces[b]->get_dfds(sTrial, aTrial));
    R.Append(temp, 0);
    // ----- R[6]
    R[6]=aTrial-aConvg;  // -kappa(fSurfaces[b]->get_dfds(sTrial, aTrial), dg);
    // ----- R[7- ]
    for (unsigned i = 0;i < activeSurfaces.size();i++)
      R[7+i]=-(fSurfaces[activeSurfaces[i]]->get_f(sTrial, aTrial));

    //=====================================================================
    // Step 4: Check convergence
    //=====================================================================
    // for (unsigned i = 0;i < fSurfaces.size();i++)
    //  cout << k<<'\t'<<i<<'\t'<<fSurfaces[i]->get_f(sTrial, aTrial)<<endl;
    // cout << R.Twonorm()<<std::endl;
    bool converged = true;
    if (R.Twonorm()>tol1) converged = false;
    for (unsigned i = 0;i < activeSurfaces.size();i++)
      if (fSurfaces[activeSurfaces[i]]->get_f(sTrial, aTrial)>tol2)
        converged = false;
    if (converged) break;

    //=====================================================================
    // Step 5: Compute elastic moduli and consistent tangent moduli
    //=====================================================================
    // Matrix A
    Matrix A(7+activeSurfaces.size(), 7+activeSurfaces.size(), 0.);
    // ----- A[0, 0-5, 5]
    Matrix A1(6, 6, 0.);
    A1.Append(invCel, 0, 0);
    for (unsigned i = 0;i < activeSurfaces.size();i++)
        A1+=dg[i]*(fSurfaces[activeSurfaces[i]]->get_df2dss(sTrial, aTrial));
    A.Append(A1, 0, 0);
    // ----- A[0, 6-6, 6]
    Vector V1(6, 0.);
    for (unsigned i = 0;i < activeSurfaces.size();i++)
        V1+=dg[i]*(fSurfaces[activeSurfaces[i]]->get_df2dsa(sTrial, aTrial));
    A.AppendCol(V1, 0, 6);
    // ----- A[0, 7-0, n]
    for (unsigned i = 0; i < activeSurfaces.size(); i++) {
      V1 = fSurfaces[activeSurfaces[i]]->get_dfds(sTrial, aTrial);
      A.AppendCol(V1, 0, 7+i);
    }
    // A[6, 0-6, 5]
    for (int i = 0;i < 6;i++) A(6, i)=0.;
    // A[6, 6-6, 6]
    A(6, 6)=1.;
    // A[6, 7] - A[7, 7+n]
    for (unsigned i = 0;i < activeSurfaces.size();i++) A(6, 7+i)=0.;
    // A[7, 0-7+n, 7+n]
    for (unsigned i = 0; i < activeSurfaces.size(); i++) {
      V1 = fSurfaces[activeSurfaces[i]]->get_dfds(sTrial, aTrial);
      A.AppendRow(V1, 7+i, 0);
    }
    // cout << sTrial << std::endl;
    // cout << sTrial.eigenvalues()<<std::endl;
    // A.report();
    // cout << "Det : "<<det(A)<<std::endl;;

    Vector x(7+activeSurfaces.size(), 0.);
    A.Solve(x, R);

    // =========================================================================
    // Step 6: Obtain increment to consistency parameter
    // =========================================================================
    Vector dEp(6, 0.);
    Vector ddg(12, 0.);
    // Stresses
    for (int i = 0;i < 6;i++) dEp[i]=x[i];
    dEp=-invCel*dEp;
    // Hardening variables
    double da = x[6];

    bool reset = false;
    // counter = 0;
    for (unsigned i = 0;i < activeSurfaces.size();i++)
      ddg[i]=x[7+i];
    for (unsigned i = 0; i < activeSurfaces.size(); i++) {
      if (dg[i]+ddg[i] < 0.) {
        fSurfaces[activeSurfaces[i]]->set_active(false);
        dg[i]=0.;
        ddg[i]=0.;
        reset = true;
      }
    }
    if (reset) continue;

    // =========================================================================
    // Step 8: Update
    // =========================================================================
    for (int i = 0;i < 6;i++) sTrial[i]+=x[i];
    ePTrial+=dEp;
    aTrial+=da;
    for (unsigned i = 0;i < activeSurfaces.size();i++)
      dg[i]+=ddg[i];
  }
  if (k == nIter) {
    // cout << "FAILED : "<<activeSurfaces.size()<<std::endl;
    // cout << sTrial << std::endl;
    // cout << ss << std::endl;
    aTrial  = aConvg;
    ePTrial = ePConvg;
    eTrial  = eTotal+De;
    sTrial  = Cel*(eTrial-ePTrial);
    Vector dg(12, 0.);
  }
}
