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

#include "domain/domain.h"

/**
 * Default Constructor.
 */
Domain::Domain()   {
  dim_ = 0;
  myTag = TAG_DOMAIN_NOT_SET;
  myFac = 1.0;
  upToDate = false;
  theDatabase = 0;
  this->init();
}
/**
 * Destructor.
 */
Domain::~Domain()  {
  this->clear();
  Containers::map_delete(theGroups);
  if (theDatabase != 0) delete theDatabase;
}
void Domain::init() {
  // Define group 0
  Group* pGroup = new Group(0);
  pGroup->set_domain(this);
  this->add(this->get_groups(), pGroup);
  // Define that groups are by set by materials
  groupsByMaterial = true;
  currentGroup = 0;
  // Initialize time/lambda
  timeCurr = 0.;
  timePrev = 0.;
  lambdaConvg = 0.;
  // Gravity direction vector and default orientation/acceleration
  gravityVect.resize(3, 0.);
  gravityVect[1]=-1;
  gravityAccl = 9.81;
}
/**
 * Clear the Domain.
 * The containers and the contained objects are deleted.
 */
void Domain::clear() {
  dim_ = 0;
  RayleighFactors.resize(0);
  eigenVals.resize(0);
  Containers::map_delete(theNodes);
  Containers::map_delete(theElements);
  Containers::map_delete(theGroups);
  Containers::map_delete(theCrossSections);
  Containers::map_delete(theMaterials);
  Containers::map_delete(theLoadCases);
  Containers::map_delete(theConstraints);
  this->init();
}
/**
 * Set the number of the dimensions of the Domain.
 * @return An integer implying if something is wrong. 
 */
int Domain::set_dim(int dim) {
  if (dim_ != 0) {
    throw SException("[nemesis:%d] %s", 9999,
      "Domain must be cleared before changing dim.");
  } else if ((dim < 1) || (dim>3)) {
    throw SException("[nemesis:%d] %s", 9999,
      "Dim value 1, 2 or 3 is only allowed.");
  } else {
    dim_ = dim;
  }
  return 0;
}
/**
 * Access to the number of the dimensions of the Domain.
 * @return The number of the dimensions of the Domain.
 */
int Domain::get_dim() const {
  return dim_;
}
void Domain::set_database(Database* pDB) {
  if (theDatabase != 0) delete theDatabase;
  theDatabase = pDB;
}
Database* Domain::get_database() {
  return theDatabase;
}
void Domain::zeroNodalStress() {
  for (NodeIterator nIter = theNodes.begin();
    nIter != theNodes.end(); nIter++)
      nIter->second->zeroStress();
}
void Domain::state(double facD) {
  for (NodeIterator nIter = theNodes.begin();
    nIter != theNodes.end(); nIter++)
      nIter->second->multDisp(facD);
}
void Domain::zeroSensitivityParameters() {
  for (ElementIterator eIter = theElements.begin();
    eIter != theElements.end(); eIter++)
      eIter->second->activateParameter(0);
}
int Domain::storeState(const char* tableName) {
  // Begin transaction
  theDatabase->beginTransaction();
  theDatabase->createTable(tableName);
  theDatabase->useTable(tableName);
  // Store nodes
  for (NodeIterator nIter = theNodes.begin();
    nIter != theNodes.end(); nIter++)
    if (nIter->second->isActive())
      theDatabase->storeData(nIter->second->get_packet());
  // Store elements
  for (ElementIterator eIter = theElements.begin();
    eIter != theElements.end(); eIter++)
    if (eIter->second->isActive())
      theDatabase->storeData(eIter->second->get_packet());
  // Commit transaction
  theDatabase->commitTransaction();
  return 0;
}
int Domain::restoreState(const char* tableName) {
  // Begin transaction
  theDatabase->beginTransaction();
  theDatabase->useTable(tableName);
  // Restore the nodes
  for (NodeIterator nIter = theNodes.begin(); nIter != theNodes.end(); nIter++)
      nIter->second->set_packet(theDatabase->retrieveData(nIter->second->get_tag(),
        nIter->second->get_id()));
  theDatabase->commitTransaction();
  return 0;
}
/**
 * Get the state of the domain.
 * A c-style string with json format is returned.
 * @return c-style string
 */
const char* Domain::get_state() {
  // define an output string stream
  ostringstream s;
  // save data
  this->save(s);
  // convert to c style string and return
  // needs to be converted to a static string before
  ///@todo: check for refactoring
  static string tmp;
  tmp=s.str();
  return tmp.c_str();
}
/**
 * Serialize domain to an output stream.
 * Follows the json format.
 */
void Domain::save(std::ostream& s) {
  s << "{";
  // Store domain info

  // Store nodes
  s << "\"nodes\":[";
  for (NodeIterator n = theNodes.begin(); n != theNodes.end(); n++) {
      if (n->second->isActive()) {
        if (n != theNodes.begin()) s << ",";
        n->second->save(s);
      }
  }
  s << "],";

  // Store elements
  s << "\"elements\":[";
  for (ElementIterator e = theElements.begin(); e != theElements.end(); e++) {  
    if (e->second->isActive()) {
      if (e != theElements.begin()) s << ",";
      e->second->save(s);
    }
  }
  s << "]";
  
  // domain
  s << "}";
}

// Rayleigh damping
void Domain::set_Rayleigh_factors(const Vector& factors) {
  RayleighFactors.resize(factors.size());
  RayleighFactors = factors;
}
const Vector& Domain::get_rayleigh_factors() {
  return RayleighFactors;
}
// EigenValues
void Domain::set_eigenvalues(const Vector& vals) {
  eigenVals.resize(vals.size());
  eigenVals = vals;
}
const Vector& Domain::get_eigen_values() {
  return eigenVals;
}
void Domain::applyLoads(double lambda_, double time_) {
  for (LoadCaseIterator lcIter = theLoadCases.begin();
    lcIter != theLoadCases.end(); lcIter++)
      lcIter->second->applyLoads(lambda_, time_);
}
void Domain::zeroLoads() {
  for (NodeIterator nIter = theNodes.begin();
    nIter != theNodes.end(); nIter++)
      nIter->second->zeroLoad();
  for (ElementIterator eIter = theElements.begin();
    eIter != theElements.end(); eIter++)
      eIter->second->zeroLoad();
}
void Domain::zeroGroups() {
  GroupIterator gi;
  for (gi = theGroups.begin(); gi != theGroups.end(); gi++)
    gi->second->set_default();
}
/**
 * Sets gravity info.
 * gravityVect Vector giving the direction of gravity. It does not include 
 * special cases for 1D, 2D or 3D cases. It might not be normalized; it is
 * normalized here and remains so thereafter.
 * @param g  The gravity acceleration.
 * @param xG x-coordinate of gravity vector.
 * @param yG y-coordinate of gravity vector.
 * @param zG z-coordinate of gravity vector.
 */
void Domain::set_gravity(double g, double xG, double yG, double zG) {
  gravityAccl = g;
  gravityVect[0]=xG;
  gravityVect[1]=yG;
  gravityVect[2]=zG;
  gravityVect.normalize();
}
/**
 * Returns gravity vector.
 * @return A constant reference to a 3x1 normalized Vector holding gravity
 * directions.
 */
const Vector& Domain::get_gravity_vect() {
  return gravityVect;
}
/**
 * Returns gravity acceleration.
 * @return Gravity acceleration
 */
double Domain::get_gravity_accl() {
  return gravityAccl;
}
