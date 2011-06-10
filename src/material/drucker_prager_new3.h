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

#ifndef SRC_MATERIAL_DRUCKER_PRAGER_NEW3_H_
#define SRC_MATERIAL_DRUCKER_PRAGER_NEW3_H_

#include <vector>
#include "material/hardening.h"
#include "material/multiaxial_material.h"
#include "numeric/matrix.h"

class YS;

/**
 * The Drucker-Prager Class.
 */
class DruckerPragerNew3: public MultiaxialMaterial {
 public:
  DruckerPragerNew3();
  DruckerPragerNew3(int id,
                   MultiaxialMaterial* elastic,
                   double c,
                   double phi,
                   double psi,
                   double Kci,
                   double Kphi,
                   double T);
  ~DruckerPragerNew3();

  MultiaxialMaterial* get_clone();
  void set_strain(const Vector& De);
  void commit();
  const Matrix& get_C();
  bool isPlastic();

  void Save(std::ostream* os);

 private:
  MultiaxialMaterial* myElastic;
  static Matrix C;
  static Matrix C3;
  bool plastic;
  int inaccurate;
  double aTrial, aConvg;
  std::vector<YS*> fSurfaces;
  std::vector<YS*> gSurfaces;
  Hardening EL;
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  DruckerPragerNew3(const DruckerPragerNew3&);
  void operator=(const DruckerPragerNew3&);
};
#endif  // SRC_MATERIAL_DRUCKER_PRAGER_NEW3_H_
