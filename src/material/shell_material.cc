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

#include "material/shell_material.h"

Matrix ShellMaterial::C(8, 8);

ShellMaterial::ShellMaterial() {
}

ShellMaterial::ShellMaterial(int id, double E, double nu, double thickness,
               double rho, double aT)
:Material(id, rho, aT) {
  sTrial.Resize(8, 0.);
  sConvg.Resize(8, 0.);
  eTrial.Resize(8, 0.);
  eTotal.Resize(8, 0.);
  
  // Material parameters
  MatParams[0] = E;
  MatParams[1] = nu;
  MatParams[2] = thickness;
}

ShellMaterial::~ShellMaterial() {
}

ShellMaterial* ShellMaterial::get_clone() {
  // Material parameters
  double E         = MatParams[ 0];
  double nu        = MatParams[ 1];
  double thickness = MatParams[ 2];
  double rho       = MatParams[30];
  double aT        = MatParams[31];
  // Create clone and return
  ShellMaterial* clone =
    new ShellMaterial(id_, E, nu, thickness, rho, aT);
  return clone;
}

void ShellMaterial::set_strain(const Vector& De) {
  eTrial = eTotal + De;
  sTrial = sConvg + (this->get_C()) * De;
}

const Matrix& ShellMaterial::get_C() {
  double E         = MatParams[ 0];
  double nu        = MatParams[ 1];
  double thickness = MatParams[ 2];
  double M  = thickness * E / (1.0 - nu * nu);
  double G  = thickness * 0.5 * E / (1.0 + nu);
  double D  =  E * (thickness * thickness * thickness) / 12.0 / ( 1.0 - nu * nu);
  C(0, 0) = M;
  C(1, 1) = M;
  C(2, 2) = G;
  C(0, 1) = nu * M;
  C(1, 0) = C(0, 1);
  G = 5. / 6. * G;
  C(3, 3) = -D;
  C(4, 4) = -D;
  C(3, 4) = -nu * D;
  C(4, 3) = C(3, 4);
  C(5, 5) = -0.5 * D * ( 1.0 - nu);
  C(6, 6) = G;
  C(7, 7) = G;
  return C;
}

void ShellMaterial::set_stress(const Vector& s) {
  sTrial = s;
  sConvg = s;
}

void ShellMaterial::addStress(const Vector& s) {
  sTrial += s;
  sConvg += s;
}

const Vector& ShellMaterial::get_stress() {
  return sTrial;
}

bool ShellMaterial::isPlastic() {
  return false;
}

void ShellMaterial::Commit() {
  eTotal = eTrial;
  sConvg = sTrial;
}

void ShellMaterial::Save(std::ostream* os) {
  // start saving
  (*os) << "{";
  (*os) << "\"data\":{";
  (*os) << "\"sigm\":"    << sConvg << ',';
  (*os) << "\"epst\":"    << eTotal << ',';
//  (*os) << "\"epsp\":"    << ePConvg << ',';
  (*os) << "\"epsv\":"    << eTotal[0]+eTotal[1]+eTotal[2] << ',';
  (*os) << "\"p\":"       << sConvg.p() << ',';
  (*os) << "\"q\":"       << sConvg.q();
  (*os) << "}";
  // finalize
  (*os) << "}";
}
