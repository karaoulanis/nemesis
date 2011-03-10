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
#include <map>

class Element;

class InitialStresses: public InitialCondition {
 public:
  InitialStresses();
  InitialStresses(const std::map<int, Element*>* elements,
                  int group_id, int dir,
                  double h1, double s1, double h2, double s2, double K0);
  ~InitialStresses();
  inline int get_group_id()   {return group_id_;}
  inline int get_dir()        {return dir_;}
  inline double get_H1()      {return h1_;}
  inline double get_S1()      {return s1_;}
  inline double get_H2()      {return h2_;}
  inline double get_S2()      {return s2_;}
  inline double get_K0()      {return K0_;}
  int Apply();

 private:
  const std::map<int, Element*>* elements_;
  int group_id_;
  int dir_;
  double h1_;
  double s1_;
  double h2_;
  double s2_;
  double K0_;
};
#endif  // SRC_LOADCASE_INITIAL_STRESSES_H_
