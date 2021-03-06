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

#ifndef SRC_MATERIAL_MATERIAL_H_
#define SRC_MATERIAL_MATERIAL_H_

#include "domain/domain_object.h"
#include "numeric/vector.h"

class Domain;

/**
 * The Material Class.
 */
class Material: public DomainObject {
 public:
  Material();
  Material(int ID, double rho, double aT);
  virtual ~Material();

  void set_X(double x1_,
             double x2_ = 0,
             double x3_ = 0);
  inline void   set_param(int i, double d)  {MatParams[i]=d;}
  inline double get_param(int i)            {return MatParams[i];}
  inline double get_rho()                   {return MatParams[30];}
  inline double get_aT()                    {return MatParams[31];}
  virtual void Commit()=0;

 protected:
  static int counter;
  int index;
  Vector MatParams;
  double x, y, z;

 private:
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  Material(const Material&);
  void operator=(const Material&);
};

#endif  // SRC_MATERIAL_MATERIAL_H_
