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

#ifndef SRC_ELEMENTS_TRIANGLE6_H_
#define SRC_ELEMENTS_TRIANGLE6_H_

#include <vector>
#include "elements/element.h"

class MatPoint;
class MultiaxialMaterial;

class Triangle6: public Element {
 public:
  // Constructors and Destructor
  Triangle6();
  Triangle6(int id,
        std::vector<Node*> nodes,
        MultiaxialMaterial* material,
        double thickness);
  ~Triangle6();

  const Matrix& get_K();
  const Matrix& get_M();
  const Vector& get_R();

  void Update(const double Dt=0.);
  void Commit();

  void get_shape_functions(MatPoint* pMatPoint, Matrix* N, double* detJ);
  void AddInitialStresses(int direction, double h1, double s1,
                          double h2, double s2, double K0);
  void recoverStresses();
  int get_num_plastic_points();

 private:
  double thickness_;
  std::vector<MatPoint*> myMatPoints;
};

#endif  // SRC_ELEMENTS_TRIANGLE6_H_
