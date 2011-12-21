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
 public:
  explicit Analysis(Domain* pDomain);
  ~Analysis();

  inline Model* get_model()                      {return &M;}
  inline AnalysisType* get_analysis_type()       {return analysis_type_;}
  inline Algorithm* get_algorithm()              {return algorithm_;}
  inline Control* get_control()                  {return control_;}
  inline Imposer* get_imposer()                  {return imposer_;}
  inline ConvergenceNorm* get_convergence_norm() {return norm_;}
  inline Reorderer* get_reorderer()              {return reorderer_;}
  inline SOE* get_soe()                          {return soe_;}
  inline Domain* get_domain()                    {return domain_;}

  void set_analysis_type(AnalysisType* analysis_type);
  void set_algorithm(Algorithm* algorithm);
  void set_control(Control* control);
  void set_imposer(Imposer* imposer);
  void set_reorderer(Reorderer* reorderer);
  void set_soe(SOE* soe);

  int Analyze(LoadCase* loadcase, int num_loadsteps);
  void Clear();

 private:
  Model M;
  Domain* domain_;
  AnalysisType* analysis_type_;
  Algorithm* algorithm_;
  Control* control_;
  Imposer* imposer_;
  ConvergenceNorm* norm_;
  Reorderer* reorderer_;
  SOE* soe_;
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  Analysis(const Analysis&);
  void operator=(const Analysis&);
};
#endif  // SRC_ANALYSIS_ANALYSIS_H_
