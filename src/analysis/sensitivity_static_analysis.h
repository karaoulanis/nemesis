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

#ifndef SRC_ANALYSIS_SENSITIVITY_STATIC_ANALYSIS_H_
#define SRC_ANALYSIS_SENSITIVITY_STATIC_ANALYSIS_H_

#include "analysis/analysis_type.h"

class Analysis;
class SensitivityControl;

class SensitivityStaticAnalysis :public AnalysisType {
 public:
  SensitivityStaticAnalysis();
  ~SensitivityStaticAnalysis();
  bool CheckIfAllows(FEObject* f);
  int Run(LoadCase* loadcase, int num_loadsteps);

 private:
  SensitivityControl* theSensitivityControl;
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  SensitivityStaticAnalysis(const SensitivityStaticAnalysis&);
  void operator=(const SensitivityStaticAnalysis&);
};

#endif  // SRC_ANALYSIS_SENSITIVITY_STATIC_ANALYSIS_H_
