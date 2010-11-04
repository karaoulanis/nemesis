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

#ifndef SRC_PARSER_PARSER_H_
#define SRC_PARSER_PARSER_H_

#include "algorithm/bfgs.h"
#include "algorithm/linear_algorithm.h"
#include "algorithm/newton_raphson_full.h"
#include "algorithm/newton_raphson_initial.h"
#include "algorithm/newton_raphson_modified.h"
#include "analysis/analysis.h"
#include "analysis/eigen_analysis.h"
#include "analysis/sensitivity_static_analysis.h"
#include "analysis/static_analysis.h"
#include "analysis/transient_analysis.h"
#include "analysis/xfem_analysis.h"
#include "constraints/constraint.h"
#include "control/arc_length_spherical.h"
#include "control/arc_length_unp.h"
#include "control/displacement_control.h"
#include "control/load_control.h"
#include "control/newmark.h"
#include "crack/crack.h"
#include "domain/domain.h"
#include "elements/bar2s.h"
#include "elements/bar2t.h"
#include "elements/beam2e.h"
#include "elements/brick8b.h"
#include "elements/brick8d.h"
#include "elements/brick8i.h"
#include "elements/quad4b.h"
#include "elements/quad4d.h"
#include "elements/quad4_disp_plain.h"
#include "elements/quad4_disp_axisymmetric.h"
#include "elements/quad4e.h"
#include "elements/quad4i.h"
#include "elements/sdof_element.h"
#include "elements/spring.h"
#include "elements/tetrahedron4_disp.h"
#include "elements/timoshenko2d.h"
#include "elements/triangle3.h"
#include "elements/triangle3_xfem.h"
#include "elements/triangle6.h"
#include "group/group.h"
#include "imposer/elimination_imposer.h"
#include "imposer/penalty_imposer.h"
#include "imposer/lagrange_imposer.h"
#include "loadcase/element_sensitivity_parameter.h"
#include "loadcase/ground_motion_file.h"
#include "loadcase/ground_motion_sin.h"
#include "loadcase/initial_displacement.h"
#include "loadcase/initial_velocity.h"
#include "loadcase/initial_stresses.h"
#include "loadcase/nodal_load_constant.h"
#include "loadcase/nodal_load_linear.h"
#include "loadcase/nodal_load_sin.h"
#include "loadcase/uniaxial_load.h"
#include "material/creep.h"
#include "material/drucker_prager.h"
#include "material/drucker_prager_new.h"
#include "material/drucker_prager_new2.h"
#include "material/drucker_prager_new3.h"
#include "material/duncan_chang.h"
#include "material/hoek_brown.h"
#include "material/lade_duncan.h"
#include "material/modified_cam_clay.h"
#include "material/mohr_coulomb.h"
#include "material/multiaxial_elastic.h"
#include "material/plane_stress.h"
#include "material/sdof_material.h"
#include "material/spring_elastic.h"
#include "material/spring_contact.h"
#include "material/spring_material.h"
#include "material/tresca.h"
#include "material/uniaxial_cyclic.h"
#include "material/uniaxial_elastic.h"
#include "material/uniaxial_elastic_plastic.h"
#include "material/uniaxial_gap.h"
#include "material/von_mises.h"
#include "model/model.h"
#include "node/node.h"
#include "reorderer/forward_cuthill_mckee.h"
#include "reorderer/forward_sloan.h"
//#include "reorderer/king.h"
//#include "reorderer/minimum_degree_ordering.h"
#include "reorderer/reverse_cuthill_mckee.h"
#include "reorderer/reverse_sloan.h"
#include "reorderer/reorderer.h"
#include "soe/full_linear_soe.h"
#include "soe/symm_linear_soe.h"
#include "soe/band_linear_soe.h"

class Parser {
 protected:
  Domain D;
  Analysis A;    
  public:
  // Constructors and destructor
  Parser();
  virtual ~Parser();

  // Parse the problem
  virtual int parse()=0;
  virtual int parse(char* filename)=0;
};

#endif  // SRC_PARSER_PARSER_H_
