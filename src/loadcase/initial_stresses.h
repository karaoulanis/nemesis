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

#ifndef SRC_LOADCASE_INITIAL_STRESSES_H_
#define SRC_LOADCASE_INITIAL_STRESSES_H_

#include "loadcase/initial_condition.h"
#include "elements/element.h"

typedef std::map < int, Element*>  ElementContainer;

class InitialStresses: public InitialCondition {
  private:
  int theGroupID;
  int dir;
  double h1;
  double s1;
  double h2;
  double s2;
  double K0;
  public:
  InitialStresses();
  InitialStresses(int groupID_, int dir_, double h1_, double s1_, double h2_, double s2_, double K0_);
  ~InitialStresses();
  inline int getGroupID() {return theGroupID;}
  inline int getDir()   {return dir;}
  inline double getH1() {return h1;}
  inline double getS1() {return s1;}
  inline double getH2() {return h2;}
  inline double getS2() {return s2;}
  inline double getK0() {return K0;}
  int apply();
};

#endif  // SRC_LOADCASE_INITIAL_STRESSES_H_
