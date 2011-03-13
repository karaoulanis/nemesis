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

#ifndef SRC_ELEMENTS_QUAD4_H_
#define SRC_ELEMENTS_QUAD4_H_

#include <vector>
#include "elements/element.h"

class MatPoint;

class Quad4: public Element {
 protected:
  int p1;                 // Integration rule xi
  int p2;                 // Integration rule eta
  static Matrix N;
  static double detJ;
  std::vector<MatPoint*> myMatPoints;
  public:
  // Constructors and Destructor
  Quad4();
  Quad4(int ID, int Node_1, int Node_2, int Node_3, int Node_4, int matID,
         int integrationRuleXi, int integrationRuleEta);
  ~Quad4();

  void commit();

  void findShapeFunctionsAt(MatPoint* pMatPoint);
  bool checkIfAllows(FEObject* f);
  void AddInitialStresses(int direction, double h1, double s1,
                          double h2, double s2, double K0);
  void recoverStresses();

  // Tracker member functions
  void addTracker(int index);
  Tracker* get_tracker(int index);
  void track();

  int get_num_plastic_points();
};

#endif  // SRC_ELEMENTS_QUAD4_H_
