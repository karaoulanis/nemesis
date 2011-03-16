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


#include "model/penalty_model_element.h"
#include "constraints/constraint.h"

/**
 * Default constructor
 */
PenaltyModelElement::PenaltyModelElement()
  :ModelElement() {
}
/**
 * Constructor
 */
PenaltyModelElement::PenaltyModelElement(const IDContainer& FTable,
                       Constraint* pConstraint,
                       double aFactor)
  :ModelElement(FTable, 0, pConstraint) {
  a = aFactor;
  myMatrix = theStaticMatrices[FTable.size()];
  myVector = theStaticVectors[FTable.size()];
}
/**
 * Destructor
 */
PenaltyModelElement::~PenaltyModelElement() {
}
void PenaltyModelElement::add_K(double factor) {
  for (unsigned int i = 0; i < theFTable.size(); i++) {
    for (unsigned int j = 0; j < theFTable.size(); j++) {
      double ci = myConstraint->get_cdof(i).coeff;
      double cj = myConstraint->get_cdof(j).coeff;
      (*myMatrix)(i, j)+=factor*a*ci*cj;
    }
  }
}
void PenaltyModelElement::add_M(double factor) {
  for (unsigned int i = 0; i < theFTable.size(); i++) {
    for (unsigned int j = 0; j < theFTable.size(); j++) {
      double ci = myConstraint->get_cdof(i).coeff;
      double cj = myConstraint->get_cdof(j).coeff;
      (*myMatrix)(i, j)+=factor*a*ci*cj;
    }
  }
}
void PenaltyModelElement::add_C(double /*factor*/) {
}
void PenaltyModelElement::add_R(double factor) {
  // Find constraint violation
  double c = myConstraint->get_val(0.);
  double cv = 0.;
  for (unsigned i = 0; i < theFTable.size(); i++) {
    double ci = myConstraint->get_cdof(i).coeff;
    double ui = myConstraint->get_disp(i);
    cv+=a*ci*ui;
  }
  cv-=a*c;
  // Add to constraint load vector
  for (unsigned i = 0; i < theFTable.size(); i++) {
    double ci = myConstraint->get_cdof(i).coeff;
    (*myVector)[i]=+factor*ci*cv;
  }
}
void PenaltyModelElement::add_Reff(double factor) {
  // Find constraint violation
  double c = myConstraint->get_val(0.);
  double cv = 0.;
  for (unsigned i = 0; i < theFTable.size(); i++) {
    double ci = myConstraint->get_cdof(i).coeff;
    double ui = myConstraint->get_disp(i);
    cv+=a*ci*ui;
  }
  cv-=a*c;
  // Add to constraint load vector
  for (unsigned i = 0; i < theFTable.size(); i++) {
    double ci = myConstraint->get_cdof(i).coeff;
    (*myVector)[i]=+factor*ci*cv;
  }
}
