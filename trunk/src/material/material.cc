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

// Included files
#include "material/material.h"

// Variable initialization
int Material::counter = 0;

/**
 * Default constructor.
 */
Material::Material() {
}
/**
 * Constructor.
 * Creates a Material and increases Material's index by one.
 * @param ID Material id.
 * @param rho Density.
 * @param aT Thermal coefficient.
 */
Material::Material(int ID, double rho, double aT)
:DomainObject(ID) {
  index = counter++;
  MatParams.resize(32);
  MatParams[30]=rho;
  MatParams[31]=aT;
}

/**
 * Destructor.
 */
Material::~Material() {
}

/**
 * Set global coordinates to material.
 * This is useful in case where material is position depended.
 */
void Material::set_X(double x1_, double x2_, double x3_) {
  x = x1_;
  y = x2_;
  z = x3_;
}
