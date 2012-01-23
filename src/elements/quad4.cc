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

#include "elements/quad4.h"
#include "material/matpoint.h"
#include "material/multiaxial_material.h"
#include "node/node.h"

Matrix Quad4::N(3, 4);
double Quad4::detJ;

Quad4::Quad4()
    : p1(0),
      p2(0),
      thickness_(0),
      myMatPoints(0),
      materials(0) {
}

Quad4::Quad4(int id, std::vector<Node*> nodes, MultiaxialMaterial* material,
             double thickness)
    : Element(id, nodes),
      p1(2),
      p2(2),
      thickness_(thickness),
      myMatPoints(4),
      materials(0) {

  // Get nodal data
  myNodalIDs.resize(4);
  myNodalIDs[0] = nodes_[0]->get_id();
  myNodalIDs[1] = nodes_[1]->get_id();
  myNodalIDs[2] = nodes_[2]->get_id();
  myNodalIDs[3] = nodes_[3]->get_id();

  // Set local nodal dofs
  myLocalNodalDofs.resize(2);
  myLocalNodalDofs[0] = 0;
  myLocalNodalDofs[1] = 1;

  // Handle common info: Start -------------------------------------------------
  // Find own matrix and vector
  myMatrix = theStaticMatrices[8];
  myVector = theStaticVectors[8];
  // Get nodal coordinates
  /// @todo replace with const references
  x.Resize(4, 3);
  for (int i = 0; i < 4; i++) {
    x(i, 0) = nodes_[i]->get_x1();
    x(i, 1) = nodes_[i]->get_x2();
    x(i, 2) = nodes_[i]->get_x3();
  }
  // Inform the nodes that the corresponding Dof's must be activated
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 2; j++)
      nodes_[i]->addDofToNode(myLocalNodalDofs[j]);
  // Load vector
  P.Resize(8, 0.);
  // Self weight
  G.Resize(8, 0.);
  myMaterial = material;
  this->AssignGravityLoads();
  // Handle common info: End ---------------------------------------------------

  // Matpoints
  myMatPoints[0]=new MatPoint(material, 1, 1, p1, p2);
  myMatPoints[1]=new MatPoint(material, 2, 1, p1, p2);
  myMatPoints[2]=new MatPoint(material, 2, 2, p1, p2);
  myMatPoints[3]=new MatPoint(material, 1, 2, p1, p2);
  for (unsigned i = 0; i < 4; i++) {
    this->findShapeFunctionsAt(myMatPoints[i]);
    double xG = N(0, 0)*x(0, 0)+N(0, 1)*x(1, 0)+N(0, 2)*x(2, 0)+N(0, 3)*x(3, 0);
    double yG = N(0, 0)*x(0, 1)+N(0, 1)*x(1, 1)+N(0, 2)*x(2, 1)+N(0, 3)*x(3, 1);
    myMatPoints[i]->set_X(xG, yG);
  }
}

Quad4::~Quad4() {
  Containers::vector_delete(&myMatPoints);
  Containers::vector_delete(&materials);
}

void Quad4::findShapeFunctionsAt(MatPoint* pMatPoint) {
  double xi= pMatPoint->get_r();
  double eta = pMatPoint->get_s();

  N(0, 0)=0.25*(1-xi)*(1-eta);           // N1
  N(0, 1)=0.25*(1+xi)*(1-eta);           // N2
  N(0, 2)=0.25*(1+xi)*(1+eta);           // N3
  N(0, 3)=0.25*(1-xi)*(1+eta);           // N4

  Matrix J(2, 2);
  J(0, 0) = 0.25*(-x(0, 0)*(1-eta)+x(1, 0)*(1-eta)
                  +x(2, 0)*(1+eta)-x(3, 0)*(1+eta));
  J(0, 1) = 0.25*(-x(0, 0)*(1-xi) -x(1, 0)*(1+xi)
                  +x(2, 0)*(1+xi) +x(3, 0)*(1-xi));
  J(1, 0) = 0.25*(-x(0, 1)*(1-eta)+x(1, 1)*(1-eta)
                  +x(2, 1)*(1+eta)-x(3, 1)*(1+eta));
  J(1, 1) = 0.25*(-x(0, 1)*(1-xi) -x(1, 1)*(1+xi)
                  +x(2, 1)*(1+xi) +x(3, 1)*(1-xi));

  detJ = J(0, 0)*J(1, 1)-J(0, 1)*J(1, 0);     // detJ
  if (detJ < 0.) {
    throw SException("[nemesis:%d] Element %d determinant negative",
                     9999, this->get_id());
  }

  double dxidx  = J(1, 1)/detJ;
  double detadx =-J(1, 0)/detJ;
  double dxidy  =-J(0, 1)/detJ;
  double detady = J(0, 0)/detJ;

  N(1, 0)=-0.25*(1-eta)*dxidx -0.25*(1-xi)*detadx;  // N1, 1
  N(1, 1)=+0.25*(1-eta)*dxidx -0.25*(1+xi)*detadx;  // N2, 1
  N(1, 2)=+0.25*(1+eta)*dxidx +0.25*(1+xi)*detadx;  // N3, 1
  N(1, 3)=-0.25*(1+eta)*dxidx +0.25*(1-xi)*detadx;  // N4, 1

  N(2, 0)=-0.25*(1-eta)*dxidy -0.25*(1-xi)*detady;  // N1, 2
  N(2, 1)=+0.25*(1-eta)*dxidy -0.25*(1+xi)*detady;  // N2, 2
  N(2, 2)=+0.25*(1+eta)*dxidy +0.25*(1+xi)*detady;  // N3, 2
  N(2, 3)=-0.25*(1+eta)*dxidy +0.25*(1-xi)*detady;  // N4, 2
}
void Quad4::Commit() {
  for (unsigned int i = 0; i < myMatPoints.size(); i++)
    myMatPoints[i]->get_material()->Commit();
}


void Quad4::AddInitialStresses(int direction,
                                 double h1, double s1,
                                 double h2, double s2, double K0) {
  for (unsigned i = 0; i < myMatPoints.size(); i++) {
    myMatPoints[i]->AddInitialStresses(direction, h1, s1, h2, s2, K0);
  }
}

void Quad4::recoverStresses() {
  /// @todo check
  if (p1 != 2||p2 != 2) return;
  double sq3 = 1.7320508075688772935274463415059;
  static Vector sigma(6);
  static Matrix E(4, 4);
  E(0, 0) =  1.0+0.5*sq3;
  E(0, 1) = -0.5;
  E(0, 2) =  1.0-0.5*sq3;
  E(0, 3) = -0.5;
  E(1, 0) = -0.5;
  E(1, 1) =  1.0+0.5*sq3;
  E(1, 2) = -0.5;
  E(1, 3) =  1.0-0.5*sq3;
  E(2, 0) =  1.0-0.5*sq3;
  E(2, 1) = -0.5;
  E(2, 2) =  1.0+0.5*sq3;
  E(2, 3) = -0.5;
  E(3, 0) = -0.5;
  E(3, 1) =  1.0-0.5*sq3;
  E(3, 2) = -0.5;
  E(3, 3) =  1.0+0.5*sq3;
  sigma.Clear();
  for (unsigned i = 0; i < 4; i++) {
    sigma[0]=E(i, 0)*myMatPoints[0]->get_material()->get_stress()[0]
        +E(i, 1)*myMatPoints[1]->get_material()->get_stress()[0]
        +E(i, 2)*myMatPoints[2]->get_material()->get_stress()[0]
        +E(i, 3)*myMatPoints[3]->get_material()->get_stress()[0];
    sigma[1]=E(i, 0)*myMatPoints[0]->get_material()->get_stress()[1]
        +E(i, 1)*myMatPoints[1]->get_material()->get_stress()[1]
        +E(i, 2)*myMatPoints[2]->get_material()->get_stress()[1]
        +E(i, 3)*myMatPoints[3]->get_material()->get_stress()[1];
    sigma[2]=E(i, 0)*myMatPoints[0]->get_material()->get_stress()[2]
        +E(i, 1)*myMatPoints[1]->get_material()->get_stress()[2]
        +E(i, 2)*myMatPoints[2]->get_material()->get_stress()[2]
        +E(i, 3)*myMatPoints[3]->get_material()->get_stress()[2];
    sigma[3]=E(i, 0)*myMatPoints[0]->get_material()->get_stress()[3]
        +E(i, 1)*myMatPoints[1]->get_material()->get_stress()[3]
        +E(i, 2)*myMatPoints[2]->get_material()->get_stress()[3]
        +E(i, 3)*myMatPoints[3]->get_material()->get_stress()[3];
    nodes_[i]->addStress(sigma);
  }
}

int Quad4::get_num_plastic_points() {
  int n = 0;
  for (unsigned int i = 0; i < myMatPoints.size(); i++) {
    if (myMatPoints[i]->get_material()->isPlastic()) n++;
  }
  return n;
}
