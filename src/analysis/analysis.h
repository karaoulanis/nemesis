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

#ifndef SRC_ANALYSIS_ANALYSIS_H_
#define SRC_ANALYSIS_ANALYSIS_H_

#include "algorithm/algorithm.h"
#include "analysis/analysis_type.h"
#include "control/control.h"
#include "convergence/convergence_norm.h"
#include "domain/domain.h"
#include "imposer/imposer.h"
#include "model/model.h"
#include "reorderer/reorderer.h"
#include "soe/soe.h"

class Algorithm;
class AnalysisType;
class Control;
class ConvergenceNorm;
class Imposer;
class Model;
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
  Analysis(Domain* pDomain);
  ~Analysis();

  inline Model* getModel()                      {return &M;}
  inline AnalysisType* getAnalysisType()        {return theAnalysisType;}
  inline Algorithm* getAlgorithm()              {return theAlgorithm;}
  inline Control* getControl()                  {return theControl;}
  inline Imposer* getImposer()                  {return theImposer;}
  inline ConvergenceNorm* getConvergenceNorm()  {return theNorm;}
  inline Reorderer* getReorderer()              {return theReorderer;}
  inline SOE* getSOE()                          {return theSOE;}
  inline Domain* getDomain()                    {return theDomain;}

  void setAnalysisType(AnalysisType* p);
  void setAlgorithm(Algorithm* p);
  void setControl(Control* p);
  void setImposer(Imposer* p);
  void setReorderer(Reorderer* p);
  void setSOE(SOE* p);

  int analyze(int lcID, int nLoadSteps);
  void clear();
};
#endif  // SRC_ANALYSIS_ANALYSIS_H_
