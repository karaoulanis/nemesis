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

#ifndef SRC_ELEMENTS_SPRING_H_
#define SRC_ELEMENTS_SPRING_H_

#include "elements/element.h"
#include "numeric/matrix.h"

class SpringMaterial;

class Spring: public Element {
 public:
  // Constructors and Destructor
  Spring();
  Spring(int id,
         std::vector<Node*> nodes,
         SpringMaterial* material,
         int dim,
         double xp1 = 1.,
         double xp2 = 0.,
         double xp3 = 0.,
         double yp1 = 0.,
         double yp2 = 1.,
         double yp3 = 0.);
  ~Spring();
  void update();
  void commit();
  bool checkIfAllows(FEObject* f);
  void recoverStresses();

  const Matrix& get_M();
  const Matrix& get_K();
  const Vector& get_R();
  const Vector& get_Reff();
  const Vector& get_Rgrad();

  // Tracker member functions
  void addTracker(int index);
  Tracker* get_tracker(int index);
  void track();

 protected:
  int dim_;
  // double gap;
  SpringMaterial* mySpringMaterial;
  Matrix T;

 private:
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  Spring(const Spring&);
  void operator=(const Spring&);
};
#endif  // SRC_ELEMENTS_SPRING_H_
