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

#ifndef SRC_MATERIAL_MULTIAXIAL_MATERIAL_H_
#define SRC_MATERIAL_MULTIAXIAL_MATERIAL_H_

#include "material/material.h"
#include "numeric/matrix.h"
#include "numeric/vector.h"

/**
 * The Multiaxial Material Class.
 */
class MultiaxialMaterial: public Material {
 protected:
  Vector sTrial;
  Vector sConvg;
  Vector eTrial;
  Vector eTotal;
  static Matrix C;
 public:
  MultiaxialMaterial();
  MultiaxialMaterial(int ID, double rho, double aT);
  ~MultiaxialMaterial();

  virtual MultiaxialMaterial* get_clone()=0;
  virtual void set_strain(const Vector& De, const double Dt=0.)=0;
  virtual const Matrix& get_C()=0;
  void set_stress(const Vector& s);
  void addStress(const Vector& s);
  const Vector& get_stress();
  virtual bool isPlastic();
};
#endif  // SRC_MATERIAL_MULTIAXIAL_MATERIAL_H_
