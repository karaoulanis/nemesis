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

#include "material/spring_material.h"
#include <sstream>
#include <string>

SpringMaterial::SpringMaterial() {
}
SpringMaterial::SpringMaterial(int ID)
:Material(ID, 0., 0.) {
  // Initialize vectors
  nDim = pD->get_dim();
  sTrial.resize(nDim, 0.);
  sConvg.resize(nDim, 0.);
  eTrial.resize(nDim, 0.);
  eTotal.resize(nDim, 0.);
  Ct.resize(3, 3, 0.);
}
/**
 */
const Matrix& SpringMaterial::get_C() {
  return Ct;
}
/**
 */
void SpringMaterial::commit() {
  sConvg = sTrial;
  eTotal = eTrial;
  this->track();
}
/**
 * Add a record to the tracker.
 * If \a myTracker pointer is null (no tracker is added) just return.
 * Otherwise gather info and send them to the tracker.
 * The domain should be already updated!
 * @todo
 */
void SpringMaterial::track() {
  if (myTracker == 0) return;
  // define an output string stream
  std::ostringstream s;
  // start saving
  s << "{";
  // save lambda
  s << "\"lambda\":"  << pD->get_lambda() << ",";
  // save time
  s << "\"time\":"    << pD->get_time_curr() << ",";
  // save self
  // save self
  s << "\"data\":{";
  s << "\"sigm\":"    << sConvg << ',';
  s << "\"epst\":"    << eTotal << ',';
//  s << "\"epsp\":"    << ePConvg << ',';
  s << "\"epsv\":"    << eTotal[0]+eTotal[1]+eTotal[2] << ',';
  s << "\"p\":"       << sConvg.p() << ',';
  s << "\"q\":"       << sConvg.q();
  s << "}";
  // finalize
  s << "}";
  // convert to c style string and return
  // needs to be converted to a static string before
  /// @todo: check for refactoring
  static string tmp;
  tmp = s.str();
  myTracker->track(tmp.c_str());
}
