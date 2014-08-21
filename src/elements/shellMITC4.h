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

#ifndef SRC_ELEMENTS_SHELLMITC4_H_
#define SRC_ELEMENTS_SHELLMITC4_H_

#include <vector>
#include "elements/element.h"

class MatPoint;
class ShellMaterial;

class ShellMITC4: public Element {
 public:
  // Constructors and Destructor
  ShellMITC4();
  ShellMITC4(int id,
        std::vector<Node*> nodes,
        ShellMaterial* material);
  ~ShellMITC4();

  const Matrix& get_K();
  const Matrix& get_M();
  const Vector& get_R();

  void Update(const double Dt=0.);
  void Commit();

  void AddInitialStresses(int direction, double h1, double s1,
                          double h2, double s2, double K0);
  void recoverStresses();

  void FindShapeFunctions();
  void FindBasis();

  void FormB(Matrix* B, const Matrix& Bm, const Matrix& Bb, const Matrix& Bs);
  void FormBm(Matrix* Bm, int node, int gp);
  void FormBb(Matrix* Bb, int node, int gp);
  void FormBd(Vector* Bd, int node, int gp);
  void FormBs(Matrix* Bs, int node, int gp);

 protected:
  std::vector<ShellMaterial*> materials_;
  static double shp_[3][4][4];
  static double detJ_[4];
  double gps_[4][2];
  double Ktt_;
  Vector v1_;
  Vector v2_;
  Vector v3_;
  double xl_[2][4];
};
#endif  // SRC_ELEMENTS_SHELLMITC4_H_
