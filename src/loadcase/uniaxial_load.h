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

#ifndef SRC_LOADCASE_UNIAXIAL_LOAD_H_
#define SRC_LOADCASE_UNIAXIAL_LOAD_H_

#include "loadcase/elemental_load.h"

/******************************************************************************
* Uniaxial Load
******************************************************************************/
class UniaxialLoad: public ElementalLoad {
 protected:
  Vector myDirection;
  Vector projections;
  double L;
  public:
  UniaxialLoad();
  UniaxialLoad(int elemID, const char* dir);
  ~UniaxialLoad();
};

/******************************************************************************
* Beam load point
******************************************************************************/
class BeamLoadPoint: public UniaxialLoad {
 private:
  double a0;
  double p0;
 public:
  BeamLoadPoint();
  BeamLoadPoint(int elemID, const char* dir, double a0, double p0);
  ~BeamLoadPoint();
  const Vector& get_P();
};

/******************************************************************************
* Beam load uniform
******************************************************************************/
class BeamLoadUniform: public UniaxialLoad {
 private:
  double p0;
 public:
  BeamLoadUniform();
  BeamLoadUniform(int elemID, const char* dir, double p0);
  ~BeamLoadUniform();
  const Vector& get_P();
};
#endif  // SRC_LOADCASE_UNIAXIAL_LOAD_H_
