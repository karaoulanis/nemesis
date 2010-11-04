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

#ifndef SRC_MATERIAL_MATERIAL_H_
#define SRC_MATERIAL_MATERIAL_H_

#include "domain/domain.h"
#include "domain/domain_object.h"
#include "numeric/vector.h"
#include "tracker/tracker.h"

class Domain;

/**
 * The Material Class.                                                
 */
class Material: public DomainObject {
 protected:
  static int counter;
  int index;
  Vector MatParams;
  Tracker* myTracker;
  double x, y, z;
  public:
  Material();
  Material(int ID, double rho, double aT);
  virtual ~Material();

  void setX(double x1_, double x2_ = 0, double x3_ = 0);
  inline void   setParam(int i, double d)  {MatParams[i]=d;}
  inline double getParam(int i)     {return MatParams[i];}
  inline double getRho()          {return MatParams[30];}
  inline double getaT()         {return MatParams[31];}
  virtual void commit()=0;

  // Tracker member functions
  void addTracker();
  Tracker* getTracker();
  virtual void track();
};

#endif  // SRC_MATERIAL_MATERIAL_H_
