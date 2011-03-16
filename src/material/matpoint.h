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

#ifndef SRC_MATERIAL_MATPOINT_H_
#define SRC_MATERIAL_MATPOINT_H_

#include "domain/domain_object.h"

class MultiaxialMaterial;

class MatPoint: public DomainObject {
  private:
  double x, y, z;
  double r, s, t;
  double w;
  static const double GaussCrds[7][7];
  static const double GaussWght[7][7];
  MultiaxialMaterial* myMaterial;
  static int IDCounter;
  public:
  MatPoint();
  MatPoint(MultiaxialMaterial* mat, int index, int p1);
  MatPoint(MultiaxialMaterial* mat, int index1, int index2, int p1, int p2);
  MatPoint(MultiaxialMaterial* mat, int index1, int index2, int index3,
                                                    int p1, int p2, int p3);
  MatPoint(MultiaxialMaterial* mat, double r_, double s_, double t_, double w_);
  ~MatPoint();

  void set_X(double x1_, double x2_ = 0, double x3_ = 0);
  inline MultiaxialMaterial* get_material()  {return myMaterial;}
  bool isPlastic();
  void AddInitialStresses(int direction, double h1, double s1,
                          double h2, double s2, double K0);
  inline double get_x()     {return x;}
  inline double get_y()     {return y;}
  inline double get_z()     {return z;}
  inline double get_r()     {return r;}
  inline double get_s()     {return s;}
  inline double get_t()     {return t;}
  inline double get_w()     {return w;}
  const Packet& get_packet()   {return thePacket;  }
  void set_packet(const Packet& /*p*/) {/*does nothing*/ }
};

#endif  // SRC_MATERIAL_MATPOINT_H_
