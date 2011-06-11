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

#include "loadcase/element_sensitivity_parameter.h"
#include "elements/element.h"

int ElementSensitivityParameter::nSensitivityParameters = 0;

ElementSensitivityParameter::ElementSensitivityParameter()
    : DomainObject(),
      element_(0),
      parameter_(0) {
}


ElementSensitivityParameter::ElementSensitivityParameter(Element* element,
                                                         int parameter)
    : DomainObject(++nSensitivityParameters),
      element_(element),
      parameter_(parameter) {
}


int ElementSensitivityParameter::apply() {
  element_->activateParameter(parameter_);
  return 0;
}

void ElementSensitivityParameter::Save(std::ostream* /*os*/) {
  /// @todo Implement this method.
}

