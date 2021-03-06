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

#ifndef SRC_MATERIAL_MOHR_COULOMB_H_
#define SRC_MATERIAL_MOHR_COULOMB_H_

#include "material/multiaxial_elastic_plastic.h"

/**
 * The Mohr-Coulomb Class.
 */
class MohrCoulomb: public MultiaxialMaterial {
 public:
  MohrCoulomb();
  MohrCoulomb(int id,
              MultiaxialMaterial* elastic,
              double c,
              double phi,
              double alpha,
              double eta);
  ~MohrCoulomb();

  MultiaxialMaterial* get_clone();
  void set_strain(const Vector& De, const double Dt=0.);
  void Commit();
  const Matrix& get_C();
  bool isPlastic();

  void Save(std::ostream* os);

 private:
  MultiaxialMaterial* myElastic;
  static Matrix C;
  static Matrix C3;
  bool plastic;
  int inaccurate;
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  MohrCoulomb(const MohrCoulomb&);
  void operator=(const MohrCoulomb&);
};
#endif  // SRC_MATERIAL_MOHR_COULOMB_H_
