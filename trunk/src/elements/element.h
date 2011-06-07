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

#ifndef SRC_ELEMENTS_ELEMENT_H_
#define SRC_ELEMENTS_ELEMENT_H_

// Included Files
#include <vector>
#include "containers/containers.h"
#include "domain/domain_object.h"
#include "numeric/matrix.h"
#include "numeric/vector.h"

// Forward Declerations
class Material;
class Node;
class Tracker;
struct GroupData;

/**
 * The Element Class.
 * The Element class is an abstract class, that provides the interface to all
 * kind of elements. One, two or three dimensional elements should be able to
 * be derived from this class.
 */
class Element: public DomainObject {
 public:
  /**
   * Default constructor.
   */
  Element();

  /**
   * Constructor.
   */
  Element(int ID, int matID);

  /**
   * Constructor.
   * @param id Element id.
   * @param nodes Pointers to element nodes.
   */
  Element(int id, std::vector<Node*> nodes);

  /**
   * Destructor.
   */
  ~Element();

  const IDContainer& get_nodal_ids() const;
  const IDContainer& get_local_nodal_dofs() const;
  const std::vector<Node*>& get_nodes() const;

  virtual int get_num_plastic_points();

  // Build and return global element matrices
  virtual const Matrix& get_K()=0;
  virtual const Matrix& get_M()=0;
  virtual const Matrix& get_C();
  virtual const Vector& get_R()=0;
  virtual const Vector& get_Reff();
  virtual const Vector& get_Rgrad();

  // Handle elemental loads
  void addLoad(const Vector& val, double fac = 1.0);
  void zeroLoad();
  void addGroundMotion(int dof, double val);
  virtual void AddInitialStresses(int direction, double h1, double s1,
                                  double h2, double s2, double K0);

  virtual void update()=0;
  virtual void commit()=0;

  virtual const Vector& get_disp_trial();
  virtual const Vector& get_velc_trial();
  virtual const Vector& get_accl_trial();
  virtual const Vector& get_disp_convg();
  virtual const Vector& get_velc_convg();
  virtual const Vector& get_accl_convg();
  virtual const Vector& get_disp_incrm();

  // GroupData related info
  void SetGroupData(const GroupData* groupdata);
  bool IsActive();
  virtual void recoverStresses() {}
  // Serialise
  void save(std::ostream& s);
  void activateParameter(int param) {activeParameter = param;}
  // Tracker member functions
  virtual void addTracker(int index);
  virtual Tracker* get_tracker(int index);
  virtual void track();

  // Enrichment functions
  virtual void enrich();

  // Temporary functions till pointers to domain are removed
  const Vector& get_gravity_vect();
  const double  get_gravity_accl();

 protected:
  std::vector<Node*> nodes_;
  IDContainer myNodalIDs;
  IDContainer myLocalNodalDofs;
  const GroupData* groupdata_;

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

 private:
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  Element(const Element&);
  void operator=(const Element&);
};
#endif  // SRC_ELEMENTS_ELEMENT_H_
