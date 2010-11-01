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

#ifndef NEMESIS_SOE_FULL_EIGEN_SOE_H_
#define NEMESIS_SOE_FULL_EIGEN_SOE_H_

#include "soe/soe.h"

class EigenSOE: public SOE {
 protected:
  valarray < double > M;
  valarray < double > ALPHAR;
  valarray < double > ALPHAI;
  valarray < double > BETA;
  valarray < double > VL;
  valarray < double > VR;
  valarray < double > WORK;
  public:
  // Constructor and destructor
  EigenSOE();
  ~EigenSOE();

  int insertMatrixIntoA(const Matrix& Ke, const IDContainer& EFTable, 
              double factor = 1.0);
  int insertMatrixIntoM(const Matrix& Ke, const IDContainer& EFTable, 
              double factor = 1.0);
  
  void zeroM();
  int solve();
  void setTheSize();
  void print();
  int getEigenSign();
};
#endif  // NEMESIS_SOE_FULL_EIGEN_SOE_H_
