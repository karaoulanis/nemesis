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

#ifndef SRC_MATERIAL_MULTIAXIAL_ELASTIC_PLASTIC_H_
#define SRC_MATERIAL_MULTIAXIAL_ELASTIC_PLASTIC_H_

#include <vector>
#include "material/multiaxial_material.h"
#include "numeric/matrix.h"
#include "numeric/vector.h"

class Surface;
class EvolutionLaw;

/**
 * The MultiaxialElastoPlastic Class.
 */
class MultiaxialElastoPlastic: public MultiaxialMaterial {
 public:
  MultiaxialElastoPlastic();
  MultiaxialElastoPlastic(int id, MultiaxialMaterial* elastic);
  ~MultiaxialElastoPlastic();

  void set_strain(const Vector& De, const double Dt=0.);
  void Commit();
  const Matrix& get_C();
  bool isPlastic()              {return plastic;}

  void updateStateVariable();

  void Save(std::ostream* os);

 protected:
  MultiaxialMaterial* myElastic;

  static Matrix C;
  Vector ePTrial, ePConvg;
  Vector qTrial, qConvg;
  double aTrial, aConvg;
  bool plastic;
  int nHardeningVariables;

  std::vector<Surface*> fSurfaces;
  std::vector<Surface*> gSurfaces;
  EvolutionLaw* EL;

  inline std::vector<Surface*> get_fSurfaces()   {return fSurfaces;}
  inline std::vector<Surface*> get_gSurfaces()   {return gSurfaces;}
  void returnMapSYS(const Vector& De);
  void returnMapSYS2(const Vector& De);

  void returnMapTest(const Vector& De);
  void returnMapMYS(const Vector& De);
  void returnMapMYS2(const Vector& De);
  void returnMapMYS3(const Vector& De);

 private:
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  MultiaxialElastoPlastic(const MultiaxialElastoPlastic&);
  void operator=(const MultiaxialElastoPlastic&);
};
#endif  // SRC_MATERIAL_MULTIAXIAL_ELASTIC_PLASTIC_H_
