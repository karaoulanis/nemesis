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

#ifndef SRC_MATERIAL_SHELL_MATERIAL_H_
#define SRC_MATERIAL_SHELL_MATERIAL_H_

#include "material/material.h"
#include "numeric/matrix.h"
#include "numeric/vector.h"

/**
 * The Shell Material Class.
 */
class ShellMaterial: public Material {
 public:
  ShellMaterial();
  ShellMaterial(int id, double E, double nu, double thickness, double rho,
                double aT);
  ~ShellMaterial();

  ShellMaterial* get_clone();
  void set_strain(const Vector& De);
  const Matrix& get_C();
  void set_stress(const Vector& s);
  void addStress(const Vector& s);
  const Vector& get_stress();
  bool isPlastic();
  void Commit();
  void Save(std::ostream* os);
  
 private:
  Vector sTrial;
  Vector sConvg;
  Vector eTrial;
  Vector eTotal;
  static Matrix C;
};
#endif  // SRC_MATERIAL_SHELL_MATERIAL_H_
