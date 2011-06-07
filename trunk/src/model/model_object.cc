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

#include "model/model_object.h"

Vector** ModelObject::theStaticVectors;
Matrix** ModelObject::theStaticMatrices;
static bool allocatedArrays = false;

/**
 * Default constructor.
 */
ModelObject::ModelObject()
    : theFTable(0),
      myMatrix(0),
      myVector(0) {
}

/**
 * Constructor.
 */
ModelObject::ModelObject(const IDContainer& FTable)
    : theFTable(FTable),
      myMatrix(0),
      myVector(0) {
  // Create static matrices and static vectors
  if (allocatedArrays == false) {
    theStaticVectors = new Vector*[64];
    theStaticMatrices = new Matrix*[64];
    for (int i = 0; i < 64; i++) {
      theStaticVectors[i]=new Vector(i, 0.);
      theStaticMatrices[i]=new Matrix(i, i, 0.);
    }
    allocatedArrays = true;
  }
}
/**
 * Destructor.
 */
ModelObject::~ModelObject() {
  if (allocatedArrays == true) {
    for (int i = 0; i < 64; i++) {
      delete theStaticVectors[i];
      delete theStaticMatrices[i];
    }
    delete[] theStaticVectors;
    delete[] theStaticMatrices;
    allocatedArrays = false;
  }
}
/**
 * Get the freedom table.
 * @return A reference to the FTable.
 */
const IDContainer& ModelObject::get_FTable() const {
  return theFTable;
}
void ModelObject::set_FTable(const IDContainer& FTable) {
  theFTable = FTable;
}
void ModelObject::set_FTable(int index, int val) {
  theFTable[index]=val;
}
