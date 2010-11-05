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

#ifndef SRC_MATERIAL_TRESCA_H_
#define SRC_MATERIAL_TRESCA_H_

#include "material/multiaxial_material.h"

/**
 * The Tresca Material Class.                                                
 */
class Tresca: public MultiaxialMaterial {
 protected:
  MultiaxialMaterial* myElastic;
  static Matrix C;
  static Matrix C3;
  bool plastic;
  int inaccurate;
  public:
  Tresca();
  Tresca(int ID, int elasticID, double sy, double kx, double ky, double kz);
  ~Tresca();

  MultiaxialMaterial* get_clone();
  void set_strain(const Vector& De);
  void commit();
  const Matrix& get_C();
  bool isPlastic();

  // Tracker member functions
  void track();
};
#endif  // SRC_MATERIAL_TRESCA_H_
