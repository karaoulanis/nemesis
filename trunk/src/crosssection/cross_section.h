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

#ifndef NEMESIS_CROSSSECTION_CROSS_SECTION_H_
#define NEMESIS_CROSSSECTION_CROSS_SECTION_H_

#include "domain/domain_object.h"

//=============================================================
// CrossSection Class
//=============================================================
class CrossSection: public DomainObject {
  private:

 protected:
  double A;
  double As2;
  double As3;
  double J1;
  double J2;
  double J3;
  double h2;
  double h3;
  public:
  CrossSection() {}
  CrossSection(int ID)
    :DomainObject(ID)  {}
    ~CrossSection() {}
  double getA()   {return A;}
  double getAs2()   {return As2;}
  double getAs3()   {return As3;}
  double getJ1()    {return J1;}
  double getJ2()    {return J2;}
  double getJ3()    {return J3;}
  double geth2()    {return h2;}
  double geth3()    {return h3;}
};

//=============================================================
// Rectangular CrossSection
//=============================================================
class RectangularCrossSection: public CrossSection {
  private:
  double w;
  double h;
  public:
  RectangularCrossSection(int ID, double width, double height);
  ~RectangularCrossSection() {}

};
//=============================================================
// User Defined CrossSection
//=============================================================
class UserDefinedCrossSection: public CrossSection {
  public:
  UserDefinedCrossSection(int ID,
              double A_in,
              double As2_in,
              double As3_in,
              double J1_in,
              double J2_in,
              double J3_in,
              double h2_in,
              double h3_in);
  ~UserDefinedCrossSection() {}
};

#endif  // NEMESIS_CROSSSECTION_CROSS_SECTION_H_
