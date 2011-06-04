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

#include "soe/soe.h"
#include <stdio.h>
#include "model/model.h"
#include "numeric/matrix.h"

SOE::SOE() 
    : model_(0),
      size_(0),
      A(),
      X(),
      B(),
      IPIV(),
      isLUFactored(false) {
}

SOE::SOE(Model* model) 
    : model_(model),
      size_(0),
      A(),
      X(),
      B(),
      IPIV(),
      isLUFactored(false) {
}

SOE::~SOE() {
}

int SOE::insertVectorIntoB(const Vector& Ve, const IDContainer& EFTable,
                           double factor) {
  for (unsigned i = 0; i < EFTable.size(); i++)
    if (EFTable[i] >= 0) B[EFTable[i]]+=factor*Ve[i];
  return 0;
}

const Vector& SOE::get_X() {
  return X;
}

const Vector& SOE::get_B() {
  return B;
}

void SOE::addB(const Vector& v) {
  B+=v;
}

void SOE::set_model(Model* model) {
  model_ = model;
}

void SOE::set_B(const Vector& v) {
  B = v;
}

void SOE::set_X(const Vector& v) {
  X = v;
}

void SOE::zero() {
  this->zeroA();
  this->zeroB();
  this->zeroX();
}

void SOE::zeroA() {
  for (unsigned i = 0; i < A.size(); i++) A[i] = 0;;
}

void SOE::zeroB() {
  B.clear();
}

void SOE::zeroX() {
  X.clear();
}

int SOE::insertMatrixIntoA(const Matrix& /*Be*/, const IDContainer& /*EFTable*/,
                  const IDContainer& /*SFTable*/, double /*factor*/) {
  return 0;
}
