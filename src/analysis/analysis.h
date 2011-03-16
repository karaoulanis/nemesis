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

#ifndef SRC_ANALYSIS_ANALYSIS_H_
#define SRC_ANALYSIS_ANALYSIS_H_

#include "model/model.h"

class Algorithm;
class AnalysisType;
class Control;
class ConvergenceNorm;
class Imposer;
class Model;
class LoadCase;
class Reorderer;
class SOE;

class Analysis {
 private:
  Model M;
  Domain* theDomain;
  AnalysisType* theAnalysisType;
  Algorithm* theAlgorithm;
  Control* theControl;
  SOE* theSOE;
  Imposer* theImposer;
  ConvergenceNorm* theNorm;
  Reorderer* theReorderer;

 public:
  explicit Analysis(Domain* pDomain);
  ~Analysis();

  inline Model* get_model()                      {return &M;}
  inline AnalysisType* get_analysis_type()       {return theAnalysisType;}
  inline Algorithm* get_algorithm()              {return theAlgorithm;}
  inline Control* get_control()                  {return theControl;}
  inline Imposer* get_imposer()                  {return theImposer;}
  inline ConvergenceNorm* get_convergence_norm() {return theNorm;}
  inline Reorderer* get_reorderer()              {return theReorderer;}
  inline SOE* get_soe()                          {return theSOE;}
  inline Domain* get_domain()                    {return theDomain;}

  void set_analysis_type(AnalysisType* p);
  void set_algorithm(Algorithm* p);
  void set_control(Control* p);
  void set_imposer(Imposer* p);
  void set_reorderer(Reorderer* p);
  void set_soe(SOE* p);

  int analyze(LoadCase* loadcase, int num_loadsteps);
  void clear();
};
#endif  // SRC_ANALYSIS_ANALYSIS_H_
