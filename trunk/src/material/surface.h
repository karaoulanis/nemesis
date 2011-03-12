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

#ifndef SRC_MATERIAL_SURFACE_H_
#define SRC_MATERIAL_SURFACE_H_

#include "numeric/matrix.h"
#include "numeric/vector.h"

class Surface {
 protected:
  static Vector myVector;
  static Matrix myMatrix;
  static Vector a1, a2, a3;
  static Matrix da2, da3;
  static double C1, C2, C3, C4;
  static double C11, C23, C32, C22, C33;
  void find_a(const Vector& s);
  void find_da(const Vector& s);
  virtual void find_C(const Vector& s, const double a);
  bool active;
  public:
  Surface();
  virtual ~Surface();

  virtual double get_f(const Vector& s, const double a)=0;
  virtual const Vector& get_dfds(const Vector& s, const double a);
  virtual double  get_dfda(const Vector& s, const double a);
  virtual const Matrix& get_df2dss(const Vector& s, const double a);
  virtual const Vector& get_df2dsa(const Vector& s, const double a);
  virtual double  get_df2daa(const Vector& s, const double a);

  void set_active(bool b)  {active = b;}
  bool isActive()     {return active;}
};
#endif  // SRC_MATERIAL_SURFACE_H_
