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

#ifndef SRC_MATERIAL_SPRING_MATERIAL_H_
#define SRC_MATERIAL_SPRING_MATERIAL_H_

#include "material/material.h"
#include "numeric/matrix.h"
#include "numeric/vector.h"

/**
 * The Single Dof Material Class.
 */
class SpringMaterial: public Material {
 protected:
  Vector sTrial;
  Vector sConvg;
  Vector eTrial;
  Vector eTotal;
  Matrix Ct;
  int nDim;
 public:
  SpringMaterial();
  explicit SpringMaterial(int ID);

  // Get clone
  virtual SpringMaterial* get_clone()=0;
  virtual void set_strain(const Vector& De)=0;

  const Matrix& get_C();
  void commit();
  inline void set_stress(const Vector& s)  {sTrial = s;    }
  inline void addStress(const Vector& s)  {sTrial+=s;   }
  inline const Vector& get_stress()    {return sTrial; }

  // Tracker member functions
  void track();
};

#endif  // SRC_MATERIAL_SPRING_MATERIAL_H_
