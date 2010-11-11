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

#ifndef SRC_MATERIAL_MULTIAXIAL_ELASTIC_PLASTIC_H_
#define SRC_MATERIAL_MULTIAXIAL_ELASTIC_PLASTIC_H_

#include "material/multiaxial_material.h"
#include "material/linear_equivalent_el.h"
#include "material/surface.h"
/**
 * The MultiaxialElastoPlastic Class.
 */
class MultiaxialElastoPlastic: public MultiaxialMaterial {
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
 public:
  MultiaxialElastoPlastic();
  MultiaxialElastoPlastic(int ID, int elasticID);
  ~MultiaxialElastoPlastic();

  void set_strain(const Vector& De);
  void commit();
  const Matrix& get_C();
  bool isPlastic()              {return plastic;}

  void updateStateVariable();

  // Tracker member functions
  void track();
};
#endif  // SRC_MATERIAL_MULTIAXIAL_ELASTIC_PLASTIC_H_
