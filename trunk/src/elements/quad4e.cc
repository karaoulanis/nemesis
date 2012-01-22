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

#include "elements/quad4e.h"
#include "group/group_data.h"
#include "main/nemesis_debug.h"
#include "material/matpoint.h"
#include "material/multiaxial_material.h"

Matrix Quad4e::Bu(3, 8);
Matrix Quad4e::Be(3, 4);
Matrix Quad4e::T0(3, 3);
Vector Quad4e::Nu(2, 8);

Quad4e::Quad4e()
  : alpha(4, 0.) {
}

Quad4e::Quad4e(int id, std::vector<Node*> nodes, MultiaxialMaterial* material,
               double thickness)
    : Quad4(id, nodes, material, thickness),
      alpha(4, 0.) {
}

Quad4e::~Quad4e() {
}

const Matrix& Quad4e::get_K() {
  this->formKR();
  Matrix &K=*myMatrix;
  double facK = groupdata_->active ? groupdata_->factor_K : 1e-7;
  K*=facK;
  return K;
}
const Matrix& Quad4e::get_M() {
  Matrix &M=*myMatrix;
  M.Clear();
  return M;
}
const Vector& Quad4e::get_R() {
  Vector& R=*myVector;
  R.Clear();
  // Quick return if inactive
  if (!(groupdata_->active)) {
    return R;
  }

  // Factors
  // double facS = groupdata_->factor_S;
  // double facG = groupdata_->factor_G;
  double facP = groupdata_->factor_P;

  this->formKR();
  R-=facP*P;
  return R;
}
void Quad4e::Update() {
}
void Quad4e::formKR() {
  static Vector Du(8);
  static Vector Depsilon(6);
  static Vector Depsilon3(3);
  static Vector sigma(6);
  static Vector sigma3(3);
  static Matrix C(3, 3);
  static Matrix Kee(4, 4);
  static Matrix Kue(8, 4);
  static Matrix Kuu(8, 8);
  static Vector Fe(4);
  static Vector Fu(8);
  static Vector Dalpha(4);
  static Vector dalpha(4);

  Du = this->get_disp_incrm();

  // Local Newton scheme to enforce orthogonality
  Dalpha.Clear();
  int counter = 0;
  do {
    counter++;
    // Form enhanced tagent matrix/vector
    Kee.Clear();
    Fe.Clear();
    for (unsigned i = 0; i < myMatPoints.size(); i++) {
      double xi =myMatPoints[i]->get_r();
      double eta = myMatPoints[i]->get_s();
      double detJ = this->get_J(xi, eta);
      double dV = myMatPoints[i]->get_w()*detJ;
      const Matrix& C0 = myMatPoints[i]->get_material()->get_C();
      C(0, 0) = C0(0, 0);
      C(0, 1) = C0(0, 1);
      C(0, 2) = C0(0, 3);
      C(1, 0) = C0(1, 0);
      C(1, 1) = C0(1, 1);
      C(1, 2) = C0(1, 3);
      C(2, 0) = C0(3, 0);
      C(2, 1) = C0(3, 1);
      C(2, 2) = C0(3, 3);
      this->formBe(xi, eta);
      this->formBu(xi, eta);

      Depsilon3 = Bu*Du+Be*Dalpha;
      Depsilon.Clear();
      Depsilon[0]=Depsilon3[0];
      Depsilon[1]=Depsilon3[1];
      Depsilon[3]=Depsilon3[2];

      myMatPoints[i]->get_material()->set_strain(Depsilon);
      sigma = myMatPoints[i]->get_material()->get_stress();
      sigma3[0]=sigma[0];
      sigma3[1]=sigma[1];
      sigma3[2]=sigma[3];

      Kee+=Transpose(Be)*C*Be*dV;
      Fe -=Transpose(Be)*sigma3*dV;
    }
    Kee.Solve(dalpha, Fe);
    Dalpha+=dalpha;
    if (counter>30) {
      for (unsigned i = 0; i < myMatPoints.size(); i++) {
        double xi =myMatPoints[i]->get_r();
        double eta = myMatPoints[i]->get_s();
        double detJ = this->get_J(xi, eta);
        double dV = myMatPoints[i]->get_w()*detJ;
        const Matrix& C0 = myMatPoints[i]->get_material()->get_C();
        C(0, 0) = C0(0, 0);
        C(0, 1) = C0(0, 1);
        C(0, 2) = C0(0, 3);
        C(1, 0) = C0(1, 0);
        C(1, 1) = C0(1, 1);
        C(1, 2) = C0(1, 3);
        C(2, 0) = C0(3, 0);
        C(2, 1) = C0(3, 1);
        C(2, 2) = C0(3, 3);
        this->formBe(xi, eta);
        report(detJ, "detJ", 20, 12);
        report(dV,  "dV  ", 20, 12);
        report(Be, "Be", 20, 2);
        report(C, "C", 20, 2);
      }
      report(Kee, "Kee", 20, 2);
      report(Fe, "Fe",  20, 10);
      report(Dalpha, "Da",  20, 10);
      throw SException("[nemesis:%d] %s", 9999, "Error in quad4e.\n");
    }
  } while (Fe.Twonorm() > 1e-6);
  alpha+=Dalpha;
  // report(alpha, "alpha", true, 12, 9);

  Kuu.Clear();
  Kue.Clear();
  for (unsigned i = 0; i < myMatPoints.size(); i++) {
    double xi =myMatPoints[i]->get_r();
    double eta = myMatPoints[i]->get_s();
    this->formBe(xi, eta);
    this->formBu(xi, eta);
    double detJ = this->get_J(xi, eta);
    double dV = myMatPoints[i]->get_w()*detJ;
    const Matrix& C0 = myMatPoints[i]->get_material()->get_C();
    C(0, 0) = C0(0, 0);
    C(0, 1) = C0(0, 1);
    C(0, 2) = C0(0, 3);
    C(1, 0) = C0(1, 0);
    C(1, 1) = C0(1, 1);
    C(1, 2) = C0(1, 3);
    C(2, 0) = C0(3, 0);
    C(2, 1) = C0(3, 1);
    C(2, 2) = C0(3, 3);
    Kuu+=Transpose(Bu)*C*Bu*dV;
    Kue+=Transpose(Bu)*C*Be*dV;
  }


  Fu.Clear();
  for (unsigned i = 0; i < myMatPoints.size(); i++) {
    double xi =myMatPoints[i]->get_r();
    double eta = myMatPoints[i]->get_s();
    this->formBu(xi, eta);
    double detJ = this->get_J(xi, eta);
    double dV = myMatPoints[i]->get_w()*detJ;
    sigma = myMatPoints[i]->get_material()->get_stress();
    sigma3[0]=sigma[0];
    sigma3[1]=sigma[1];
    sigma3[2]=sigma[3];
    Fu+=Transpose(Bu)*sigma3*dV;
  }

  Matrix& K=*myMatrix;
  Vector& R=*myVector;
  K = Kuu-Kue*Inverse(Kee)*Transpose(Kue);
  R = Fu -Kue*Inverse(Kee)*Fe;
  // K.report("K", 12, 3);
  // Kuu.report("Kuu", 12, 3);
  // Kue.report("Kue", 12, 3);
  // Kee.report("Kee", 12, 3);
  // cout << "------------------------------------------------"<<endl;
}

double Quad4e::get_J(double xi, double eta) {
  static Matrix J(2, 2);
  J(0, 0) = 0.25*(-x(0, 0)*(1-eta)+x(1, 0)*(1-eta)
                  +x(2, 0)*(1+eta)-x(3, 0)*(1+eta));
  J(0, 1) = 0.25*(-x(0, 0)*(1-xi) -x(1, 0)*(1+xi)
                  +x(2, 0)*(1+xi) +x(3, 0)*(1-xi));
  J(1, 0) = 0.25*(-x(0, 1)*(1-eta)+x(1, 1)*(1-eta)
                  +x(2, 1)*(1+eta)-x(3, 1)*(1+eta));
  J(1, 1) = 0.25*(-x(0, 1)*(1-xi) -x(1, 1)*(1+xi)
                  +x(2, 1)*(1+xi)  +x(3, 1)*(1-xi));

  return J(0, 0)*J(1, 1)-J(0, 1)*J(1, 0);     // detJ
}

void Quad4e::formT0(double xi, double eta) {
  static Matrix J(2, 2);
  J(0, 0) = 0.25*(-x(0, 0)*(1-eta)+x(1, 0)*(1-eta)
                  +x(2, 0)*(1+eta)-x(3, 0)*(1+eta));
  J(0, 1) = 0.25*(-x(0, 0)*(1-xi) -x(1, 0)*(1+xi)
                  +x(2, 0)*(1+xi) +x(3, 0)*(1-xi));
  J(1, 0) = 0.25*(-x(0, 1)*(1-eta)+x(1, 1)*(1-eta)
                  +x(2, 1)*(1+eta)-x(3, 1)*(1+eta));
  J(1, 1) = 0.25*(-x(0, 1)*(1-xi) -x(1, 1)*(1+xi)
                  +x(2, 1)*(1+xi) +x(3, 1)*(1-xi));

  T0(0, 0) =   J(0, 0)*J(0, 0);
  T0(0, 1) =   J(1, 0)*J(0, 1);
  T0(0, 2) = 2*J(0, 0)*J(0, 1);
  T0(1, 0) =   J(0, 1)*J(1, 0);
  T0(1, 1) =   J(1, 1)*J(1, 1);
  T0(1, 2) = 2*J(1, 0)*J(1, 1);
  T0(2, 0) =   J(0, 0)*J(1, 0);
  T0(2, 1) =   J(0, 1)*J(1, 1);
  T0(2, 2) =   J(0, 0)*J(1, 1)+J(0, 1)*J(1, 0);
}

void Quad4e::formNu(double xi, double eta) {
  N.Clear();
  N(0, 0)=0.25*(1-xi)*(1-eta);  // N1
  N(0, 2)=0.25*(1+xi)*(1-eta);  // N2
  N(0, 4)=0.25*(1+xi)*(1+eta);  // N3
  N(0, 6)=0.25*(1-xi)*(1+eta);  // N4
  N(1, 1)=N(0, 0);
  N(1, 3)=N(0, 2);
  N(1, 5)=N(0, 4);
  N(1, 7)=N(0, 6);
}
void Quad4e::formBu(double xi, double eta) {
  static Matrix J(2, 2);
  J(0, 0) = 0.25*(-x(0, 0)*(1-eta) +x(1, 0)*(1-eta)
                  +x(2, 0)*(1+eta) -x(3, 0)*(1+eta));
  J(0, 1) = 0.25*(-x(0, 0)*(1-xi)  -x(1, 0)*(1+xi)
                  +x(2, 0)*(1+xi)  +x(3, 0)*(1-xi));
  J(1, 0) = 0.25*(-x(0, 1)*(1-eta) +x(1, 1)*(1-eta)
                  +x(2, 1)*(1+eta) -x(3, 1)*(1+eta));
  J(1, 1) = 0.25*(-x(0, 1)*(1-xi)  -x(1, 1)*(1+xi)
                  +x(2, 1)*(1+xi)  +x(3, 1)*(1-xi));

  double detJ = J(0, 0)*J(1, 1)-J(0, 1)*J(1, 0);  // detJ

  double dxidx  = J(1, 1)/detJ;
  double detadx =-J(1, 0)/detJ;
  double dxidy  =-J(0, 1)/detJ;
  double detady = J(0, 0)/detJ;

  double d1dxi=-0.25*(1-eta)*dxidx -0.25*(1-xi)*detadx;  // N1, 1
  double d2dxi=+0.25*(1-eta)*dxidx -0.25*(1+xi)*detadx;  // N2, 1
  double d3dxi=+0.25*(1+eta)*dxidx +0.25*(1+xi)*detadx;  // N3, 1
  double d4dxi=-0.25*(1+eta)*dxidx +0.25*(1-xi)*detadx;  // N4, 1

  double d1deta=-0.25*(1-eta)*dxidy -0.25*(1-xi)*detady;  // N1, 2
  double d2deta=+0.25*(1-eta)*dxidy -0.25*(1+xi)*detady;  // N2, 2
  double d3deta=+0.25*(1+eta)*dxidy +0.25*(1+xi)*detady;  // N3, 2
  double d4deta=-0.25*(1+eta)*dxidy +0.25*(1-xi)*detady;    // N4, 2

  Bu(0, 0) = d1dxi;
  Bu(0, 1) = 0.;
  Bu(0, 2) = d2dxi;
  Bu(0, 3) = 0.;
  Bu(0, 4) = d3dxi;
  Bu(0, 5) = 0.;
  Bu(0, 6) = d4dxi;
  Bu(0, 7) = 0.;
  Bu(1, 0) = 0.;
  Bu(1, 1) = d1deta;
  Bu(1, 2) = 0.;
  Bu(1, 3) = d2deta;
  Bu(1, 4) = 0.;
  Bu(1, 5) = d3deta;
  Bu(1, 6) = 0.;
  Bu(1, 7) = d4deta;
  Bu(2, 0) = d1deta;
  Bu(2, 1) = d1dxi;
  Bu(2, 2) = d2deta;
  Bu(2, 3) = d2dxi;
  Bu(2, 4) = d3deta;
  Bu(2, 5) = d3dxi;
  Bu(2, 6) = d4deta;
  Bu(2, 7) = d4dxi;
}
// void Quad4e::formBe(double xi, double eta)
// {
//   Be.clear();
//   Be(0, 0)=xi;
//   Be(1, 1)=eta;
//   Be(2, 2)=xi;
//   Be(2, 3)=eta;
//   //Be(2, 4)=xi*eta;
//
//   double detJ0 = get_J(0., 0.);
//   double detJ =get_J(xi, eta);
//   this->formT0(0., 0.);
//   Be=(detJ0/detJ)*Transpose(Inverse(T0))*Be;
//
// }
void Quad4e::formBe(double xi, double eta) {
  static Matrix J(2, 2);
  static Matrix invTranJ(2, 2);
  J(0, 0)=0.25*(-x(0, 0)+x(1, 0)+x(2, 0)-x(3, 0));
  J(0, 1)=0.25*(-x(0, 0)-x(1, 0)+x(2, 0)+x(3, 0));
  J(1, 0)=0.25*(-x(0, 1)+x(1, 1)+x(2, 1)-x(3, 1));
  J(1, 1)=0.25*(-x(0, 1)-x(1, 1)+x(2, 1)+x(3, 1));

  invTranJ = Transpose(Inverse(J));
  double detJ =get_J(xi, eta);

  double d00 = invTranJ(0, 0)*xi/detJ;
  double d10 = invTranJ(1, 0)*xi/detJ;
  double d01 = invTranJ(0, 1)*eta/detJ;
  double d11 = invTranJ(1, 1)*eta/detJ;

  Be.Clear();
  Be(0, 0) = d00;
  Be(0, 1) = 0.;
  Be(0, 2) = d01;
  Be(0, 3) = 0.;
  Be(1, 0) = 0.;
  Be(1, 1) = d10;
  Be(1, 2) = 0.;
  Be(1, 3) = d11;
  Be(2, 0) = d10;
  Be(2, 1) = d00;
  Be(2, 2) = d11;
  Be(2, 3) = d01;
}
