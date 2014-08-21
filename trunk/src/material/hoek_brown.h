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

#ifndef SRC_MATERIAL_HOEK_BROWN_H_
#define SRC_MATERIAL_HOEK_BROWN_H_

#include <vector>
#include "material/multiaxial_elastic_plastic.h"

/**
 * The Hoek-Brown Class.
 */
class HoekBrown: public MultiaxialMaterial {
 public:
  HoekBrown();
  HoekBrown(int id,
            MultiaxialMaterial* elastic,
            double si,
            double sp,
            double mb,
            double mbb,
            double alpha);
  ~HoekBrown();

  MultiaxialMaterial* get_clone();
  void set_strain(const Vector& De, const double Dt=0.);
  void Commit();
  const Matrix& get_C();
  bool isPlastic();

  void find_f(const Vector& s, double q);
  void find_dfds(const Vector& s, double q);
  void find_dgds(const Vector& s, double q);
  void find_d2gdsds(const Vector& s, double q);

  void Save(std::ostream* os);

 private:
  MultiaxialMaterial* myElastic;
  static Matrix C;
  static Matrix C3;
  bool plastic;
  int inaccurate;
  Vector f;
  std::vector<Vector> dfds;
  std::vector<Vector> dgds;
  std::vector<Matrix> d2gdsds;
  double aTrial, aConvg;
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  HoekBrown(const HoekBrown&);
  void operator=(const HoekBrown&);
};
#endif  // SRC_MATERIAL_HOEK_BROWN_H_
