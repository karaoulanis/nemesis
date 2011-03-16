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

#ifndef SRC_ALGORITHM_BFGS_H_
#define SRC_ALGORITHM_BFGS_H_

#include <vector>
#include "algorithm/algorithm.h"

class Vector;

class BFGS :public Algorithm {
 private:
  int m;
  double etaMin;
  double etaMax;
  double rTol;
  int maxIter;
  bool isLineSearchActive;
  std::vector<Vector> s;
  std::vector<Vector> y;
 public:
  explicit BFGS(int m_ = 10);
  BFGS(int m_, double etaMin_, double etaMax_, double rTol_, int maxIter_);
  ~BFGS();
  int solveStep(int n);
  void lineSearch(double s0, double s1, const Vector& du);
};
#endif  // SRC_ALGORITHM_BFGS_H_
