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

#ifndef SRC_SOE_BAND_LINEAR_SOE_H_
#define SRC_SOE_BAND_LINEAR_SOE_H_

#include "soe/soe.h"

class BandLinearSOE: public SOE {
 public:
  // Constructor and destructor
  BandLinearSOE();
  explicit BandLinearSOE(Model* model);
  ~BandLinearSOE();

  int insertMatrixIntoA(const Matrix& Ke, const IDContainer& EFTable,
              double factor = 1.0);
  int solve();
  void set_size();
  void print();
  int get_eigen_sign();
 private:
  int lowerBandwidth;
  int upperBandwidth;
  int nRows;
};
#endif  // SRC_SOE_BAND_LINEAR_SOE_H_
