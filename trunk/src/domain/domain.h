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

#ifndef SRC_DOMAIN_DOMAIN_H_
#define SRC_DOMAIN_DOMAIN_H_

#include <map>
#include <iostream>
#include <utility>
#include "numeric/vector.h"

// Forward declarations
class CrossSection;
class Element;
class Material;
class Node;
class Group;
class LoadCase;
class Constraint;
class Tracker;

// Type definitions - Containers
typedef std::map<int, Node*>                          NodeContainer;
typedef std::map<int, Element*>                       ElementContainer;
typedef std::map<int, Group*>                         GroupContainer;
typedef std::map<int, CrossSection*>                  CrossSectionContainer;
typedef std::map<int, Material*>                      MaterialContainer;
typedef std::map<int, LoadCase*>                      LoadCaseContainer;
typedef std::map<int, Constraint*>                    ConstraintContainer;
typedef std::map<int, Tracker*>                       TrackerContainer;

// Type definitions - Iterators
typedef std::map<int, Node*>::const_iterator          NodeIterator;
typedef std::map<int, Element*>::const_iterator       ElementIterator;
typedef std::map<int, Group*>::const_iterator         GroupIterator;
typedef std::map<int, CrossSection*>::const_iterator  CrossSectionIterator;
typedef std::map<int, Material*>::const_iterator      MaterialIterator;
typedef std::map<int, LoadCase*>::const_iterator      LoadCaseIterator;
typedef std::map<int, Constraint*>::const_iterator    ConstraintIterator;
typedef std::map<int, Tracker*>::const_iterator       TrackerIterator;

// Domain tag
enum DomainTag {  TAG_DOMAIN_NOT_SET      =  0,
                  TAG_DOMAIN_1D           = 10,
                  TAG_DOMAIN_2D           = 20,
                  TAG_DOMAIN_PLANE_STRESS = 21,
                  TAG_DOMAIN_PLANE_STRAIN = 22,
                  TAG_DOMAIN_AXISYMMETRIC = 23,
                  TAG_DOMAIN_3D           = 30 };
/**
 * The Domain Class.
 * The Domain is where all data concerning the nodes, elements, materials, loads
 * a.s.o. are stored. As storage schemes the stl map class is used. Domain
 * provides template methods for storing, accessing and handling these oblects.
 */
class Domain {
 private:
  int dim_;
  DomainTag myTag;
  double myFac;   // factor for plane/axisymmetric problems (e.g. 1.0 rad)
  bool upToDate;
  Vector gravityVect;
  double gravityAccl;
  bool groupsByMaterial;

  NodeContainer         nodes_;
  CrossSectionContainer sections_;
  ElementContainer      elements_;
  GroupContainer        groups_;
  MaterialContainer     materials_;
  LoadCaseContainer     loadcases_;
  ConstraintContainer   constraints_;
  TrackerContainer      trackers_;

  Vector RayleighFactors;
  Vector eigenVals;
  double timeCurr;
  double timePrev;
  double lambdaConvg;
  bool export_vtk_;

  void init();
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  Domain(const Domain&);
  void operator=(const Domain&);

 public:
  // Constructors and Destructors
  Domain();
  ~Domain();

  // Access to data members
  int set_dim(int dim);
  int get_dim() const;
  void set_export_vtk(bool b);

  // Analysis/LoadCase related members
  void Initialize();
  void zeroNodalStress();
  void zeroSensitivityParameters();
  void ApplyLoads(double lambda_, double time_);
  void zeroLoads();

  void clear();
  void state(double facD);

  // Gravity axis
  void set_gravity(double g, double xG, double yG, double zG);
  const Vector& get_gravity_vect();
  double get_gravity_accl();

  /// @todo cleanup
  double get_time_curr()      {return timeCurr;}
  double get_time_prev()      {return timePrev;}
  double get_time_incr()      {return timeCurr-timePrev;}
  void incTime(double dt)     {timeCurr+=dt;}
  void Commit();
  void set_lambda(double l)   {lambdaConvg = l;}
  double get_lambda()         {return lambdaConvg;}

  void rollback()             {}

  // serialize/deserialize
  const char* get_state();
  void Save(std::ostream* s);

  void set_tag(DomainTag t)           {myTag = t;}
  DomainTag get_tag()                 {return myTag;}
  void set_fac(double fac)            {myFac = fac;}
  double get_fac()                    {return myFac;}
  bool isUpToDate()                   {return upToDate;}
  void set_uptodate()                 {upToDate = true;}
  bool areGroupsByMaterial()          {return groupsByMaterial;}
  void set_groups_by_material(bool b) {groupsByMaterial = b;}

  // Rayleigh damping
  void set_Rayleigh_factors(const Vector& factors);
  const Vector& get_rayleigh_factors();

  // EigenValues
  void set_eigenvalues(const Vector& vals);
  const Vector& get_eigen_values();
  
  // export
  void ExportVTK();

  inline NodeContainer& get_nodes()                  {return nodes_;}
  inline CrossSectionContainer& get_sections()       {return sections_;}
  inline ElementContainer& get_elements()            {return elements_;}
  inline GroupContainer& get_groups()                {return groups_;}
  inline MaterialContainer& get_materials()          {return materials_;}
  inline LoadCaseContainer& get_loadcases()          {return loadcases_;}
  inline ConstraintContainer& get_constraints()      {return constraints_;}
  inline TrackerContainer& get_trackers()            {return trackers_;}

  template<class TE, class TC> int add(TC& c, TE* e) {
    int id = e->get_id();
    std::pair<typename TC::const_iterator, bool> res =
            c.insert(std::make_pair(id, e));
    if (res.second == true) {
      upToDate = false;
    } else {
      delete e;
      throw SException("[nemesis:%d] Component %d already exists.", 9999, id);
    }
    return 0;
  }

  template<class TC> int rem(TC& c, int id) {
    typename TC::iterator p = c.find(id);
    if (p == c.end()) {
      throw SException("[nemesis:%d] Component %d not found.", 9999, id);
    } else {
      delete p->second;
      c.erase(p);
    }
    return 0;
  }

  template<class TE, class TC> TE* get(TC& c, int id) {
    typename TC::iterator p = c.find(id);
    if (p != c.end()) {
      return p->second;
    } else {
      throw SException("[nemesis:%d] Component %d does not exist.", 9999, id);
    }
  }
};

#endif  // SRC_DOMAIN_DOMAIN_H_
