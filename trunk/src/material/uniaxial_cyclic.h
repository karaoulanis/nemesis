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

#ifndef SRC_MATERIAL_UNIAXIAL_CYCLIC_H_
#define SRC_MATERIAL_UNIAXIAL_CYCLIC_H_

#include "material/uniaxial_material.h"

class UniaxialCyclic: public UniaxialMaterial {
  private:
  double sr;
  double er;
  double eTrial;
  double eConvg;
  double Et;
  bool reversed;
  public:
  UniaxialCyclic();
  UniaxialCyclic(int ID, double E, double nu, double rho, double aT,
                 double tmax, double Gmax);
  ~UniaxialCyclic();
  UniaxialMaterial* get_clone();
  void set_strain(const double De);
  double get_C();
  void commit();
};
#endif  // SRC_MATERIAL_UNIAXIAL_CYCLIC_H_
