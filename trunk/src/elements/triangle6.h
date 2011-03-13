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

#ifndef SRC_ELEMENTS_TRIANGLE6_H_
#define SRC_ELEMENTS_TRIANGLE6_H_

#include <vector>
#include "elements/element.h"

class MatPoint;

class Triangle6: public Element {
 protected:
  double a1, a2, a3;
  double b1, b2, b3;
  double c1, c2, c3;
  double A;
  std::vector<MatPoint*> myMatPoints;
  public:
  // Constructors and Destructor
  Triangle6();
  Triangle6(int ID, int Node_1, int Node_2, int Node_3,
    int Node_4, int Node_5, int Node_6, int matID);
  ~Triangle6();

  const Matrix& get_K();
  const Matrix& get_M();
  const Vector& get_R();

  void update();
  void commit();

  void get_shape_functions(MatPoint* pMatPoint, Matrix& N, double& detJ);
  void AddInitialStresses(int direction, double h1, double s1,
                          double h2, double s2, double K0);
  bool checkIfAllows(FEObject* f);
  void recoverStresses();
  int get_num_plastic_points();
};

#endif  // SRC_ELEMENTS_TRIANGLE6_H_
