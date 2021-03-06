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

#ifndef SRC_MODEL_MODEL_OBJECT_H_
#define SRC_MODEL_MODEL_OBJECT_H_

#include "elements/element.h"

class ModelObject {
 public:
  // Constructors
  ModelObject();
  explicit ModelObject(const IDContainer& FTable);
  virtual ~ModelObject();
  // Access to data members
  const IDContainer& get_FTable() const;
  void set_FTable(const IDContainer& FTable);
  void set_FTable(int index, int val);

  inline void zeroMatrix()        {myMatrix->Clear();}
  inline void zeroVector()        {myVector->Clear();}
  inline const Matrix& get_matrix() const  {return *myMatrix; }
  inline const Vector& get_vector() const  {return *myVector; }

 protected:
  IDContainer theFTable;
  static Matrix** theStaticMatrices;
  static Vector** theStaticVectors;
  Matrix* myMatrix;
  Vector* myVector;

 private:
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  ModelObject(const ModelObject&);
  void operator=(const ModelObject&);
};
#endif  // SRC_MODEL_MODEL_OBJECT_H_
