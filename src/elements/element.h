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

#ifndef SRC_ELEMENTS_ELEMENT_H_
#define SRC_ELEMENTS_ELEMENT_H_

// Included Files
#include "containers/containers.h"
#include "domain/domain.h"
#include "domain/domain_object.h"
#include "group/group.h"
#include "loadcase/elemental_load.h"
#include "loadcase/initial_stresses.h"
#include "node/node.h"
#include "numeric/matrix.h"
#include "numeric/vector.h"

// Forward Declerations
class Domain;
class ElementalLoad;
class Group;
class InitialStresses;
class Material;
class Node;
class Tracker;

/**
 * The Element Class.                                                
 * The Element class is an abstract class, that provides the interface to all
 * kind of elements. One, two or three dimensional elements should be able to 
 * be derived from this class. 
 */
class Element: public DomainObject {
 protected:
  IDContainer myNodalIDs;
  IDContainer myLocalNodalDofs;
  std::vector<Node*> myNodes;
  Group* myGroup;

  // Materials
  Material* myMaterial;

  Vector P;
  Vector G;
  Matrix x;
  Vector b;

  Matrix* myMatrix;
  Vector* myVector;
  static Matrix** theStaticMatrices;
  static Vector** theStaticVectors;
  int handleCommonInfo();
  int activeParameter;

  bool active_;

 public:
  // Constructors and Destructor
  Element();
  Element(int ID, int matID);
  ~Element();

  const IDContainer& get_nodal_ids() const;
  const IDContainer& get_local_nodal_dofs() const;
  const std::vector<Node*>& get_nodes() const;

  bool isActive();
  virtual int get_num_plastic_points();

  // Build and return global element matrices
  virtual const Matrix& get_K()=0;
  virtual const Matrix& get_M()=0;
  virtual const Matrix& get_C();
  virtual const Vector& get_R()=0;
  virtual const Vector& get_Reff()   {myVector->clear(); return *myVector;}
  virtual const Vector& get_Rgrad()  {myVector->clear(); return *myVector;}

  // Handle elemental loads
  void addLoad(const Vector& val, double fac = 1.0);
  void zeroLoad();
  void addGroundMotion(int dof, double val);
  virtual void addInitialStresses(InitialStresses* pInitialStresses);

  virtual void update()=0;
  virtual void commit()=0;

  virtual const Vector& get_disp_trial();
  virtual const Vector& get_velc_trial();
  virtual const Vector& get_accl_trial();
  virtual const Vector& get_disp_convg();
  virtual const Vector& get_velc_convg();
  virtual const Vector& get_accl_convg();
  virtual const Vector& get_disp_incrm();

  void set_group(int groupID);
  virtual void recoverStresses() {}

  // Send and receive packet
  const Packet& get_packet();
  void set_packet(const Packet& p);
  void save(std::ostream& s);

  void activateParameter(int param) {activeParameter = param;}

  // Tracker member functions
  virtual void addTracker(int index);
  virtual Tracker* get_tracker(int index);
  virtual void track();

  // Enrichment functions
  virtual void enrich();
};

#endif  // SRC_ELEMENTS_ELEMENT_H_
