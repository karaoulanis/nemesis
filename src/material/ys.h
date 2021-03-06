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

#ifndef SRC_MATERIAL_YS_H_
#define SRC_MATERIAL_YS_H_

#include <vector>

#include "main/nemesis_debug.h"
#include "numeric/matrix.h"

class YS {
 public:
  YS();
  virtual ~YS();
  void set_sigma(const Vector& s);

  // functions to be overwritten
  virtual double get_f(const Vector& sigma, const double kappa)=0;
  virtual const Vector& get_dfds(const Vector& sigma, const double kappa)=0;
  virtual const Matrix& get_d2fdsds(const Vector& sigma, const double kappa)=0;
  virtual double get_dfdk(const Vector& sigma, const double kappa)=0;
  virtual const Vector& get_f2dkds(const Vector& sigma, const double kappa)=0;
 protected:
  double s1, s2, s3;
  double I1, J2;
  Vector a, a1, a2;
  Matrix da, da2, da11, da22;
  bool active;
};
#endif  // SRC_MATERIAL_YS_H_
