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

#include "elements/beam2e.h"
#include "crosssection/cross_section.h"
#include "group/group_data.h"
#include "material/uniaxial_material.h"
#include "node/node.h"

Beam2e::Beam2e()
    : myUniMaterial(0),
      mySection(0),
      L(0.) {
}

Beam2e::Beam2e(int id, std::vector<Node*> nodes, UniaxialMaterial* material,
       CrossSection* section)
    : Element(id, nodes),
      myUniMaterial(material),
      mySection(section) {
  // Tag
  myTag = TAG_ELEM_BEAM_2D_EULER;
  // Nodal ids
  myNodalIDs.resize(2);
  myNodalIDs[0] = nodes_[0]->get_id();
  myNodalIDs[1] = nodes_[1]->get_id();
  // Local dofs
  myLocalNodalDofs.resize(3);
  myLocalNodalDofs[0] = 0;
  myLocalNodalDofs[1] = 1;
  myLocalNodalDofs[2] = 5;
  
  // Handle common info: Start -------------------------------------------------
  // Find own matrix and vector
  myMatrix = theStaticMatrices[6];
  myVector = theStaticVectors[6];
  // Get nodal coordinates
  /// @todo replace with const references
  x.Resize(2, 3);
  for (int i = 0; i < 2; i++) {
    x(i, 0) = nodes_[i]->get_x1();
    x(i, 1) = nodes_[i]->get_x2();
    x(i, 2) = nodes_[i]->get_x3();
  }
  // Inform the nodes that the corresponding Dof's must be activated
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 3; j++)
      nodes_[i]->addDofToNode(myLocalNodalDofs[j]);
  // Load vector
  P.resize(6, 0.);
  // Self weight
  G.resize(6, 0.);
  this->AssignGravityLoads();
  // Handle common info: End ---------------------------------------------------

  // Length 
  L = sqrt((x(1, 1)-x(0, 1))*(x(1, 1)-x(0, 1))
          +(x(1, 0)-x(0, 0))*(x(1, 0)-x(0, 0)));
  // Directional cosines
  cosX[0] = (x(1, 0)-x(0, 0))/L;
  cosX[1] = (x(1, 1)-x(0, 1))/L;
  // Self weight - Transform vector b to local system
  /// @todo: check this
  double A = mySection->get_A();
  double b0A = b[0]*A;
  double b1A = b[1]*A;
  b[0]= cosX[0]*b0A+cosX[1]*b1A;
  b[1]=-cosX[1]*b0A+cosX[0]*b1A;
}

/**
 * Destructor.
 */
Beam2e::~Beam2e() {
}

const Matrix& Beam2e::get_K() {
  Matrix& K=*myMatrix;
  K.Clear();
  double E =myUniMaterial->get_param(0);
  double A = mySection->get_A();
  double J = mySection->get_J3();

  double c = cosX[0];
  double s = cosX[1];
  double C1 = E*A/L;
  double C2 = 12*E*J/(L*L*L);
  double C3 = 6*E*J/(L*L);
  double K1 = C1*c*c+C2*s*s;
  double K2=(C1-C2)*c*s;
  double K3 = C1*s*s+C2*c*c;
  double K4=-C3*s;
  double K5 = C3*c;
  double K6 = 4*E*J/L;
  double K7 = 2*E*J/L;

  K(0, 0)= K1;
  K(0, 1)= K2;
  K(0, 2)= K4;
  K(0, 3)=-K1;
  K(0, 4)=-K2;
  K(0, 5)= K4;
  K(1, 0)= K2;
  K(1, 1)= K3;
  K(1, 2)= K5;
  K(1, 3)=-K2;
  K(1, 4)=-K3;
  K(1, 5)= K5;
  K(2, 0)= K4;
  K(2, 1)= K5;
  K(2, 2)= K6;
  K(2, 3)=-K4;
  K(2, 4)=-K5;
  K(2, 5)= K7;
  K(3, 0)=-K1;
  K(3, 1)=-K2;
  K(3, 2)=-K4;
  K(3, 3)= K1;
  K(3, 4)= K2;
  K(3, 5)=-K4;
  K(4, 0)=-K2;
  K(4, 1)=-K3;
  K(4, 2)=-K5;
  K(4, 3)= K2;
  K(4, 4)= K3;
  K(4, 5)=-K5;
  K(5, 0)= K4;
  K(5, 1)= K5;
  K(5, 2)= K7;
  K(5, 3)=-K4;
  K(5, 4)=-K5;
  K(5, 5)= K6;

  double facK = groupdata_->active ? groupdata_->factor_K : 1e-7;
  K*=facK;
  return K;
}

const Matrix& Beam2e::get_M() {
  Matrix& M=*myMatrix;
  M.Clear();
  return M;
}

const Vector& Beam2e::get_Rgrad() {
  /// @todo
  Matrix& K=*myMatrix;
  K.Clear();
  myVector->clear();
  double E =myUniMaterial->get_param(0);
  double A = mySection->get_A();
  double J = mySection->get_J3();

  if (activeParameter == 1) {
    E = 1;
  } else if (activeParameter == 2) {
    A = 1;
  } else if (activeParameter == 3) {
    J = 1;
  } else {
    return *myVector;
  }

  double c = cosX[0];
  double s = cosX[1];
  double C1 = E*A/L;
  double C2 = 12*E*J/(L*L*L);
  double C3 = 6*E*J/(L*L);
  double K1 = C1*c*c+C2*s*s;
  double K2=(C1-C2)*c*s;
  double K3 = C1*s*s+C2*c*c;
  double K4=-C3*s;
  double K5 = C3*c;
  double K6 = 4*E*J/L;
  double K7 = 2*E*J/L;

  K(0, 0)= K1;
  K(0, 1)= K2;
  K(0, 2)= K4;
  K(0, 3)=-K1;
  K(0, 4)=-K2;
  K(0, 5)= K4;
  K(1, 0)= K2;
  K(1, 1)= K3;
  K(1, 2)= K5;
  K(1, 3)=-K2;
  K(1, 4)=-K3;
  K(1, 5)= K5;
  K(2, 0)= K4;
  K(2, 1)= K5;
  K(2, 2)= K6;
  K(2, 3)=-K4;
  K(2, 4)=-K5;
  K(2, 5)= K7;
  K(3, 0)=-K1;
  K(3, 1)=-K2;
  K(3, 2)=-K4;
  K(3, 3)= K1;
  K(3, 4)= K2;
  K(3, 5)=-K4;
  K(4, 0)=-K2;
  K(4, 1)=-K3;
  K(4, 2)=-K5;
  K(4, 3)= K2;
  K(4, 4)= K3;
  K(4, 5)=-K5;
  K(5, 0)= K4;
  K(5, 1)= K5;
  K(5, 2)= K7;
  K(5, 3)=-K4;
  K(5, 4)=-K5;
  K(5, 5)= K6;

  double facK = groupdata_->active ? groupdata_->factor_K: 1e-7;
  K*=facK;
  return *myVector;
}

const Vector& Beam2e::get_R() {
  Vector& R=*myVector;
  R.clear();
  // Quick return if inactive
  if (!(groupdata_->active)) {
    return R;
  }
  // Factors
  double facS = groupdata_->factor_S;
  double facG = groupdata_->factor_G;
  double facP = groupdata_->factor_P;

  double c = cosX[0];
  double s = cosX[1];
  static Vector u(6);
  u = this->get_disp_incrm();

  // Global u to local u
  double u0= cosX[0]*u[0]+cosX[1]*u[1];
  double u1=-cosX[1]*u[0]+cosX[0]*u[1];
  double u3= cosX[0]*u[3]+cosX[1]*u[4];
  double u4=-cosX[1]*u[3]+cosX[0]*u[4];
  u[0]=u0;
  u[1]=u1;
  u[3]=u3;
  u[4]=u4;

  // Find factors of Ke
  double E =myUniMaterial->get_param(0);
  double A = mySection->get_A();
  double J = mySection->get_J3();
  double KN = E*A/L;
  double K2 = 2*E*J/L;
  double K4 = 4*E*J/L;
  double K6 = 6*E*J/(L*L);
  double K12 = 12*E*J/(L*L*L);

  // R = Fint - Fext(Gravity) - Fext(Loads)
  R[0] = facS*(+KN*u[0]  - KN*u[3])                      - facG*(+0.5*b[0]*L);
  R[1] = facS*(+K12*u[1] + K6*u[2] - K12*u[4] + K6*u[5]) - facG*(+0.5*b[1]*L);
  R[2] = facS*(+K6*u[1]  + K4*u[2] - K6*u[4]  + K2*u[5]) - facG*(+b[1]*L*L/12.);
  R[3] = facS*(-KN*u[0]  + KN*u[3])                      - facG*(+0.5*b[0]*L);
  R[4] = facS*(-K12*u[1] - K6*u[2] + K12*u[4] - K6*u[5]) - facG*(+0.5*b[1]*L);
  R[5] = facS*(+K6*u[1]  + K2*u[2] - K6*u[4]  + K4*u[5]) - facG*(-b[1]*L*L/12.);
  R-=facP*P;

  // Local R to global R
  double R0 = c*R[0]-s*R[1];
  double R1 = s*R[0]+c*R[1];
  double R3 = c*R[3]-s*R[4];
  double R4 = s*R[3]+c*R[4];
  R[0]=R0;
  R[1]=R1;
  R[3]=R3;
  R[4]=R4;
  return R;
}

void Beam2e::recoverStresses() {
  return;
}
