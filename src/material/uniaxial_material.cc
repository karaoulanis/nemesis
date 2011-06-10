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


#include "material/uniaxial_material.h"
#include <sstream>
#include <string>

UniaxialMaterial::UniaxialMaterial() {
}

UniaxialMaterial::UniaxialMaterial(int ID, double rho, double aT)
    : Material(ID, rho, aT),
      sTrial(0.),
      sConvg(0.),
      eTotal(0.) {
}

UniaxialMaterial::~UniaxialMaterial() {
}

void UniaxialMaterial::Save(std::ostream* os) {
  // start saving
  (*os) << "{";
  (*os) << "\"data\":{";
  (*os) << "\"sigm\":"    << sConvg << ',';
  (*os) << "\"epst\":"    << eTotal;
//  (*os) << "\"epsp\":"    << ePConvg << ',';
//  (*os) << "\"epsv\":"    << eTotal[0]+eTotal[1]+eTotal[2] << ',';
//  (*os) << "\"p\":"       << sConvg.p() << ',';
//  (*os) << "\"q\":"       << sConvg.q();
  (*os) << "}";
  // finalize
  (*os) << "}";
}
