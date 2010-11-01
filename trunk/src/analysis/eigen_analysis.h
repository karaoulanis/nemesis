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

#ifndef NEMESIS_ANALYSIS_EIGEN_ANALYIS_H_
#define NEMESIS_ANALYSIS_EIGEN_ANALYIS_H_

#include "analysis/analysis_type.h"
#include "control/eigen_control.h"
#include "imposer/elimination_imposer.h"
#include "soe/eigen_soe.h"

class Analysis;
class EigenAnalysis :public AnalysisType {
  public:
  EigenAnalysis();
  ~EigenAnalysis();
  bool checkIfAllows(FEObject* f);
  int run(int nLC, int nLoadSteps);
};

#endif  // NEMESIS_ANALYSIS_EIGEN_ANALYIS_H_
