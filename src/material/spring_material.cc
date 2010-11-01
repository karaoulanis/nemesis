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

#include "material/spring_material.h"

SpringMaterial::SpringMaterial() {
}
SpringMaterial::SpringMaterial(int ID)
:Material(ID, 0., 0.) {
  // Initialize vectors
  nDim = pD->getnDim();
  sTrial.resize(nDim, 0.);
  sConvg.resize(nDim, 0.);
  eTrial.resize(nDim, 0.);
  eTotal.resize(nDim, 0.);
  Ct.resize(3, 3, 0.);
}
/**
 */
const Matrix& SpringMaterial::getC() {
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
  ostringstream s;
  s << "DATA "  <<' ';
//  s << "sigm "  <<1020<<' '<<sConvg<<' ';
//  s << "epst "  <<1020<<' '<<eTotal<<' ';
  s << "END "<<' ';
  myTracker->track(pD->getLambda(), pD->getTimeCurr(), s.str());
}
