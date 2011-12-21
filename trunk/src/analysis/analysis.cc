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

#include "analysis/analysis.h"
#include "algorithm/algorithm.h"
#include "analysis/analysis_type.h"
#include "control/control.h"
#include "convergence/convergence_norm.h"
#include "domain/domain.h"
#include "imposer/imposer.h"
#include "reorderer/reorderer.h"
#include "soe/soe.h"

Analysis::Analysis(Domain* domain)
    : M(domain),
      domain_(domain),
      analysis_type_(NULL),
      algorithm_(NULL),
      control_(NULL),
      soe_(NULL),
      imposer_(NULL),
      norm_(new ConvergenceNorm()),
      reorderer_(NULL) {
  norm_->set_analysis(this);
}

Analysis::~Analysis() {
  delete norm_;
  this->Clear();
}

int Analysis::Analyze(LoadCase* loadcase, int num_loadsteps) {
  if (analysis_type_ == 0)
    throw SException("[nemesis:%d] %s", 9999, "No analysis type set.");
  // Run the analysis
  return analysis_type_->Run(loadcase, num_loadsteps);
}

void Analysis::Clear() {
  M.clear();
  if (analysis_type_ != 0) {
    delete analysis_type_;
    analysis_type_ = 0;
  }
  if (imposer_ != 0) {
    delete imposer_;
    imposer_ = 0;
  }
  if (control_ != 0) {
    delete control_;
    control_ = 0;
  }
  if (algorithm_ != 0) {
    delete algorithm_;
    algorithm_ = 0;
  }
  if (soe_ != 0) {
    delete soe_;
    soe_ = 0;
  }
  if (reorderer_ != 0) {
    delete reorderer_;
    reorderer_ = 0;
  }
}
void Analysis::set_analysis_type(AnalysisType* analysis_type) {
  if (analysis_type_ != 0) delete analysis_type_;
  analysis_type_ = analysis_type;
}

void Analysis::set_algorithm(Algorithm* algorithm) {
  if (algorithm_ != 0) delete algorithm_;
  algorithm_ = algorithm;
}

void Analysis::set_control(Control* control) {
  if (control_ != 0) delete control_;
  control_ = control;
}

void Analysis::set_imposer(Imposer* imposer) {
  if (imposer_ != 0) delete imposer_;
  imposer_ = imposer;
}

void Analysis::set_reorderer(Reorderer* reorderer) {
  if (reorderer_ != 0) delete reorderer_;
  reorderer_ = reorderer;
}

void Analysis::set_soe(SOE* soe) {
  if (soe_ != 0) delete soe_;
  soe_ = soe;
  soe_->set_model(&M);
}
