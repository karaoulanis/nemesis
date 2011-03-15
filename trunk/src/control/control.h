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

#ifndef SRC_CONTROL_CONTROL_H_
#define SRC_CONTROL_CONTROL_H_

#include "analysis/analysis_object.h"
#include "numeric/vector.h"

class ModelElement;
class ModelNode;

class Control: public AnalysisObject {
 protected:
  double lambdaTrial;
  double lambdaConvg;
  double DLambda;
  double dLambda;
  Vector qRef;
 public:
  Control();
  virtual ~Control();

  virtual void formTangent();
  virtual void formResidual(double factor)=0;

  // Functions to build Element by Element or Node by Node it's contribution
  virtual void formElementalTangent(ModelElement* pModelElement)=0;
  virtual void formElementalResidual(ModelElement* pModelElement,
                                     double time = 0.)=0;
  virtual void formNodalResidual(ModelNode* pModelNode)=0;

  virtual double get_lambda();
  virtual double get_time() {return 0;}  /// @todo: implement this better

  virtual void init()=0;
  virtual void predict()=0;
  virtual void correct()=0;
  virtual void commit()=0;
  virtual void rollback()=0;
};
#endif  // SRC_CONTROL_CONTROL_H_
