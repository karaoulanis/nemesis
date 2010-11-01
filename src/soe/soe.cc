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
* along with this program.  If not, see < http://www.gnu.org/licenses/>.        *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#include "soe/soe.h"

SOE::SOE() {
  theSize = 0;
  isLUFactored = false;
}
SOE::~SOE() {
}
int SOE::insertVectorIntoB(const Vector& Ve, const IDContainer& EFTable,
                           double factor) {
  for (unsigned i = 0; i < EFTable.size(); i++)
    if (EFTable[i] >= 0) B[EFTable[i]]+=factor*Ve[i];
  return 0;
}
void SOE::printSolution() {
  for (int i = 0;i < theSize; i++) cout << X[i] << endl;
}
const Vector& SOE::getX() {
  return X;
}
const Vector& SOE::getB() {
  return B;
}
void SOE::addB(const Vector& v) {
  B+=v;
}
void SOE::setB(const Vector& v) {
  B = v;
}
void SOE::setX(const Vector& v) {
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
int SOE::plotGraph(const char* s) {
  // Create the Graph
  if (theSize <= 0)
    throw SException("[nemesis:%d] %s", 1110, "Could not initialize graph.\n");
  UndirectedGraph G(theSize);
  pA->getModel()->getUndirectedGraph(G);
  // Plot the Graph
  ofstream outfile(s);
  ///@todo: the following line is commented, produces error/warnings
  // write_graphviz(outfile, G);
  return 0;
}
