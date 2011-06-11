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

#ifndef SRC_MATERIAL_UNIAXIAL_MATERIAL_H_
#define SRC_MATERIAL_UNIAXIAL_MATERIAL_H_

#include "material/material.h"

/**
 * The Uniaxial Material Class.
 */
class UniaxialMaterial: public Material {
 public:
  UniaxialMaterial();
  UniaxialMaterial(int ID, double rho, double aT);
  ~UniaxialMaterial();

  virtual UniaxialMaterial* get_clone()=0;
  virtual void set_strain(const double De)=0;
  virtual double get_C()=0;
  inline void set_stress(const double s) {sTrial = s;}
  inline void addStress(const double s)  {sTrial+= s;}  /// @todo: check
  inline double get_stress()             {return sTrial;}
  void Save(std::ostream* os);
 protected:
  double sTrial;
  double sConvg;
  double eTotal;
};

#endif  // SRC_MATERIAL_UNIAXIAL_MATERIAL_H_
