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

#ifndef SRC_MATERIAL_CREEP_H_
#define SRC_MATERIAL_CREEP_H_

#include "material/multiaxial_material.h"
/**
 * The Creep Class.
 */
class Creep: public MultiaxialMaterial {
 public:
  Creep();
  Creep(int id,
        MultiaxialMaterial* elastic,
        double A,
        double n,
        double k);
  ~Creep();

  void set_strain(const Vector& De, const double Dt=0.);
  void Commit();
  const Matrix& get_C();
  bool isPlastic()              {return false;}
  MultiaxialMaterial* get_clone();
  void Save(std::ostream* os);

 protected:
  MultiaxialMaterial* myElastic;
  static Matrix C;
  Vector eCTrial, eCConvg;
  double A;
  double n;
  double k;

 private:
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  Creep(const Creep&);
  void operator=(const Creep&);
};
#endif  // SRC_MATERIAL_CREEP_H_
