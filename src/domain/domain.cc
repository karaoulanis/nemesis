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

#include "domain/domain.h"
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "constraints/constraint.h"
#include "crosssection/cross_section.h"
#include "domain/domain_object.h"
#include "elements/element.h"
#include "exception/sexception.h"
#include "group/group.h"
#include "group/group_data.h"
#include "loadcase/loadcase.h"
#include "material/material.h"
#include "node/node.h"
#include "tracker/tracker.h"

/**
 * Default Constructor.
 */
Domain::Domain()
    : dim_(0),
      myTag(TAG_DOMAIN_NOT_SET),
      myFac(1.0),
      upToDate(false),
      gravityVect(0),
      gravityAccl(0),
      groupsByMaterial(true),
      nodes_(),
      sections_(),
      elements_(),
      groups_(),
      materials_(),
      loadcases_(),
      constraints_(),
      trackers_(),
      RayleighFactors(0),
      eigenVals(0),
      timeCurr(0.),
      timePrev(0.),
      lambdaConvg(0.) {
  this->init();
}
/**
 * Destructor.
 */
Domain::~Domain()  {
  this->clear();
  Containers::map_delete(&groups_);
}
void Domain::init() {
  // Define group 0
  Group* pGroup = new Group(0);
  this->add(this->get_groups(), pGroup);
  // Define that groups are by set by materials
  groupsByMaterial = true;
  // Initialize time/lambda
  timeCurr = 0.;
  timePrev = 0.;
  lambdaConvg = 0.;
  // Gravity direction vector and default orientation/acceleration
  gravityVect.Resize(3, 0.);
  gravityVect[1]=-1;
  gravityAccl = 9.81;
}
/**
 * Clear the Domain.
 * The containers and the contained objects are deleted.
 */
void Domain::clear() {
  dim_ = 0;
  RayleighFactors.Resize(0);
  eigenVals.Resize(0);
  Containers::map_delete(&nodes_);
  Containers::map_delete(&elements_);
  Containers::map_delete(&groups_);
  Containers::map_delete(&sections_);
  Containers::map_delete(&materials_);
  Containers::map_delete(&loadcases_);
  Containers::map_delete(&constraints_);
  Containers::map_delete(&trackers_);
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

void Domain::zeroNodalStress() {
  for (NodeIterator nIter = nodes_.begin();
    nIter != nodes_.end(); nIter++)
      nIter->second->zeroStress();
}
void Domain::state(double facD) {
  for (NodeIterator nIter = nodes_.begin();
    nIter != nodes_.end(); nIter++)
      nIter->second->multDisp(facD);
}
void Domain::zeroSensitivityParameters() {
  for (ElementIterator eIter = elements_.begin();
    eIter != elements_.end(); eIter++)
      eIter->second->activateParameter(0);
}


const char* Domain::get_state() {
  // define an output string stream
  std::ostringstream s;
  // save data
  this->Save(&s);
  // convert to c style string and return
  // needs to be converted to a static string before
  /// @todo: check for refactoring
  static string tmp;
  tmp = s.str();
  return tmp.c_str();
}


void Domain::Save(std::ostream* os) {
  (*os) << "{";
  // Store domain info

  // Store nodes
  (*os) << "\"nodes\":[";
  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
      if (n->second->IsActive()) {
        if (n != nodes_.begin()) (*os) << ",";
        n->second->Save(os);
      }
  }
  (*os) << "],";

  // Store elements
  (*os) << "\"elements\":[";
  for (ElementIterator e = elements_.begin(); e != elements_.end(); e++) {
    if (e->second->IsActive()) {
      if (e != elements_.begin()) (*os) << ",";
      e->second->Save(os);
    }
  }
  (*os) << "]";

  // domain
  (*os) << "}";
}


void Domain::Commit() {
  timePrev = timeCurr;
  for (TrackerIterator i = trackers_.begin(); i != trackers_.end(); i++) {
    i->second->Track(lambdaConvg, timeCurr);
  }
  // this->ExportVTK();
}


// Rayleigh damping
void Domain::set_Rayleigh_factors(const Vector& factors) {
  RayleighFactors.Resize(factors.get_size());
  RayleighFactors = factors;
  Element::set_rayleigh(factors);
}

const Vector& Domain::get_rayleigh_factors() {
  return RayleighFactors;
}

// EigenValues
void Domain::set_eigenvalues(const Vector& vals) {
  eigenVals.Resize(vals.get_size());
  eigenVals = vals;
}

const Vector& Domain::get_eigen_values() {
  return eigenVals;
}

void Domain::Initialize() {
  /// @todo: rename as to be consistent
  // default groupdata for all elements
  // static GroupData default_groupdata;
  // default_groupdata.active = true;
  // default_groupdata.factor_K = 1.0;
  // default_groupdata.factor_S = 1.0;
  // default_groupdata.factor_G = 1.0;
  // default_groupdata.factor_P = 1.0;

  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
    n->second->set_active(false);
  }

  for (ElementIterator e = elements_.begin(); e != elements_.end(); e++) {
//    e->second->SetGroupData(default_groupdata);
    e->second->SetActiveNodes();
  }
}

void Domain::zeroLoads() {
  for (ElementIterator e = elements_.begin(); e != elements_.end(); e++) {
    e->second->zeroLoad();
  }
  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
    n->second->zeroLoad();
  }
}

/**
 * @brief Apply loads for given \f$\lambda\f$ and \f$t\f$.
 * @param lambda Given \f$\lambda\f$.
 * @param time Given \f$t\f$.
 */
void Domain::ApplyLoads(double lambda, double time) {
  for (LoadCaseIterator iter = loadcases_.begin();
    iter != loadcases_.end(); iter++)
      iter->second->ApplyLoads(lambda, time);
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
  gravityVect[0] = xG;
  gravityVect[1] = yG;
  gravityVect[2] = zG;
  gravityVect.Normalize();
  // Check if elements are already defined
  if (elements_.size() > 0) {
    throw SException("[nemesis:%d] %s", 9999,
                     "Gravity info must me set before elements' definitions.");
  }
  // Set gravity info to elements
  Element::set_gravitydirection(gravityVect[0], gravityVect[1], gravityVect[2]);
  Element::set_gravityacceleration(g);
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

// export
void Domain::ExportVTK() {
  std::vector<std::string> data_;
  data_.push_back("# vtk DataFile Version 3.0");
  data_.push_back("Author: F.E. Karaoulanis");
  data_.push_back("ASCII");
  data_.push_back("DATASET UNSTRUCTURED_GRID");
  data_.push_back("");

  std::map<int, int> nodes_index;
  std::stringstream ss;
  
  ///@todo avoid running the same loop twice for getting the number of points
  int num_points = 0;
  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
    if (n->second->IsActive()) {
      num_points++;
    }
  }
  ss  << "POINTS " << num_points << " float" << std::endl;
  int k = 0;
  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
    if (n->second->IsActive()) {
      ss  << std::setw(12) << std::fixed << std::setprecision(6)
          << n->second->get_x1() << " "
          << std::setw(12) << std::fixed << std::setprecision(6)
          << n->second->get_x2() << " "
          << std::setw(12) << std::fixed << std::setprecision(6)
          << n->second->get_x3() << std::endl;
      nodes_index[n->second->get_id()] = k++;
    }
  }
  data_.push_back(ss.str());
  
  ///@todo avoid running the same loop twice for getting the number of points
  int num_elems = 0;
  int num_elem_nodes = 0;
  for (ElementIterator e = elements_.begin(); e != elements_.end(); e++) {
    if (e->second->IsActive()) {
      num_elems++;
      num_elem_nodes += 1 + e->second->get_nodes().size();
    }
  }

  ss.str("");
  ss  << "CELLS "
      << num_elems      << " "
      << num_elem_nodes << std::endl;
  for (ElementIterator e = elements_.begin(); e != elements_.end(); e++) {
    if (e->second->IsActive()) {
      
      if (e->second->get_nodes().size() == 8) {  // eight node bricks
        ss  << "8 "
            << std::setw(6) << nodes_index[e->second->get_nodes()[0]->get_id()] << " "
            << std::setw(6) << nodes_index[e->second->get_nodes()[1]->get_id()] << " "
            << std::setw(6) << nodes_index[e->second->get_nodes()[2]->get_id()] << " "
            << std::setw(6) << nodes_index[e->second->get_nodes()[3]->get_id()] << " "
            << std::setw(6) << nodes_index[e->second->get_nodes()[4]->get_id()] << " "
            << std::setw(6) << nodes_index[e->second->get_nodes()[5]->get_id()] << " "
            << std::setw(6) << nodes_index[e->second->get_nodes()[6]->get_id()] << " "
            << std::setw(6) << nodes_index[e->second->get_nodes()[7]->get_id()] << std::endl;
      } else if (e->second->get_nodes().size() == 4) {  // shells, quads
        ss  << "4 "
            << std::setw(6) << nodes_index[e->second->get_nodes()[0]->get_id()] << " "
            << std::setw(6) << nodes_index[e->second->get_nodes()[1]->get_id()] << " "
            << std::setw(6) << nodes_index[e->second->get_nodes()[2]->get_id()] << " "
            << std::setw(6) << nodes_index[e->second->get_nodes()[3]->get_id()] << std::endl;      
      } else {
        std::cout << "Element type not supported." << std::endl;
      }

    }
  }
  ss << std::endl;

  ss << "CELL_TYPES " << num_elems << std::endl;
  for (ElementIterator e = elements_.begin(); e != elements_.end(); e++) {
    if (e->second->IsActive()) {
      if (e->second->get_nodes().size() == 8) {  // eight node bricks
        ss << 12 << std::endl;
      } else if (e->second->get_nodes().size() == 4) {  // shells, quads
        ss <<  9 << std::endl;
      } else {
        std::cout << "Element type not supported." << std::endl;
      }
    }
  }
 
  data_.push_back(ss.str());

  ss.str("");
  ss << "CELL_DATA " << num_elems << std::endl;
  data_.push_back(ss.str());

  ss.str("");
  ss << "SCALARS material int" << std::endl;
  ss << "LOOKUP_TABLE default" << std::endl;
  for (ElementIterator e = elements_.begin(); e != elements_.end(); e++) {
    if (e->second->IsActive()) {
      ss  << e->second->get_material()->get_id() << std::endl;
    }
  }
  data_.push_back(ss.str());

  ss.str("");
  ss << "SCALARS plastic_points int" << std::endl;
  ss << "LOOKUP_TABLE default" << std::endl;
  for (ElementIterator e = elements_.begin(); e != elements_.end(); e++) {
    if (e->second->IsActive()) {
      ss  << e->second->get_num_plastic_points() << std::endl;
    }
  }
  data_.push_back(ss.str());

  ss.str("");
  ss << "POINT_DATA " << num_points << std::endl;
  data_.push_back(ss.str());

  ss.str("");
  ss << "VECTORS disp float" << std::endl;
  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
    if (n->second->IsActive()) {
      ss  << std::setw(16) << std::fixed << std::setprecision(8)
          << n->second->get_disp_convg()[0] << " "
          << std::setw(16) << std::fixed << std::setprecision(8)
          << n->second->get_disp_convg()[1] << " "
          << std::setw(16) << std::fixed << std::setprecision(8)
          << n->second->get_disp_convg()[2] << std::endl;
    }
  }
  data_.push_back(ss.str());

  ss.str("");
  ss << "SCALARS sigma_x float" << std::endl;
  ss << "LOOKUP_TABLE default" << std::endl;
  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
    if (n->second->IsActive()) {
      ss  << std::setw(16) << std::fixed << std::setprecision(8)
          << n->second->get_stress()[0]
          << std::endl;
    }
  }
  data_.push_back(ss.str());

  ss.str("");
  ss << "SCALARS sigma_y float" << std::endl;
  ss << "LOOKUP_TABLE default" << std::endl;
  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
    if (n->second->IsActive()) {
      ss  << std::setw(16) << std::fixed << std::setprecision(8)
          << n->second->get_stress()[1]
          << std::endl;
    }
  }
  data_.push_back(ss.str());

  ss.str("");
  ss << "SCALARS sigma_z float" << std::endl;
  ss << "LOOKUP_TABLE default" << std::endl;
  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
    if (n->second->IsActive()) {
      ss  << std::setw(16) << std::fixed << std::setprecision(8)
          << n->second->get_stress()[2]
          << std::endl;
    }
  }
  data_.push_back(ss.str());

  ss.str("");
  ss << "SCALARS sigma_xy float" << std::endl;
  ss << "LOOKUP_TABLE default" << std::endl;
  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
    if (n->second->IsActive()) {
      ss  << std::setw(16) << std::fixed << std::setprecision(8)
          << n->second->get_stress()[3]
          << std::endl;
    }
  }
  data_.push_back(ss.str());

  ss.str("");
  ss << "SCALARS sigma_yz float" << std::endl;
  ss << "LOOKUP_TABLE default" << std::endl;
  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
    if (n->second->IsActive()) {
      ss  << std::setw(16) << std::fixed << std::setprecision(8)
          << n->second->get_stress()[4]
          << std::endl;
    }
  }
  data_.push_back(ss.str());

  ss.str("");
  ss << "SCALARS sigma_zx float" << std::endl;
  ss << "LOOKUP_TABLE default" << std::endl;
  for (NodeIterator n = nodes_.begin(); n != nodes_.end(); n++) {
    if (n->second->IsActive()) {
      ss  << std::setw(16) << std::fixed << std::setprecision(8)
          << n->second->get_stress()[5]
          << std::endl;
    }
  }
  data_.push_back(ss.str());

  static int lc = 0;
  ss.str("");
  ss << "lc_" << std::setw(4) << std::setfill('0') << lc << ".vtk";
  std::ofstream file(ss.str().c_str());
  for (unsigned i = 0; i < data_.size(); i++) {
    file << data_[i] << std::endl;
  }
  lc++;
}
