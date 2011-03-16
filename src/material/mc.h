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

#ifndef SRC_MATERIAL_MC_H_
#define SRC_MATERIAL_MC_H_

#include "material/surface.h"
#include "numeric/vector.h"

class MC: public Surface {
 protected:
  double c;
  double phi;
  void find_C(const Vector& s, const double a);
  virtual void find_A(const Vector& s, double& A, double& dA, double& d2A)=0;
  public:
  MC();
  MC(double c_, double phi_);
  ~MC();
  double get_f(const Vector& s, const double q);
};
#endif  // SRC_MATERIAL_MC_H_
