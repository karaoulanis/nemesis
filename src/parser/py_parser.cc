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

#include "parser/py_parser.h"

// C header files
#ifdef _DEBUG
  #include <Python.h>
  #undef Py_DEBUG
#else
  #include <Python.h>
#endif
#include <stdio.h>

// C++ header files
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// additional includes
#include "algorithm/bfgs.h"
#include "algorithm/linear_algorithm.h"
#include "algorithm/newton_raphson_full.h"
#include "algorithm/newton_raphson_initial.h"
#include "algorithm/newton_raphson_modified.h"
#include "analysis/eigen_analysis.h"
#include "analysis/sensitivity_static_analysis.h"
#include "analysis/static_analysis.h"
#include "analysis/transient_analysis.h"
#include "constraints/constraint.h"
#include "control/arc_length_spherical.h"
#include "control/arc_length_unp.h"
#include "control/displacement_control.h"
#include "control/load_control.h"
#include "control/newmark.h"
#include "convergence/convergence_norm.h"
#include "crosssection/cross_section.h"
#include "database/sqlite_database.h"
#include "elements/bar2s.h"
#include "elements/bar2t.h"
#include "elements/beam2e.h"
#include "elements/brick8b.h"
#include "elements/brick8d.h"
#include "elements/brick8i.h"
#include "elements/quad4_disp_axisymmetric.h"
#include "elements/quad4_disp_plain.h"
#include "elements/quad4b.h"
#include "elements/quad4d.h"
#include "elements/quad4e.h"
#include "elements/quad4i.h"
#include "elements/sdof_element.h"
#include "elements/spring.h"
#include "elements/tetrahedron4_disp.h"
#include "elements/timoshenko2d.h"
#include "elements/triangle3.h"
#include "elements/triangle6.h"
#include "exception/sexception.h"
#include "group/group.h"
#include "imposer/elimination_imposer.h"
#include "imposer/lagrange_imposer.h"
#include "imposer/penalty_imposer.h"
#include "loadcase/element_sensitivity_parameter.h"
#include "loadcase/ground_motion_file.h"
#include "loadcase/ground_motion_sin.h"
#include "loadcase/group_state.h"
#include "loadcase/initial_displacement.h"
#include "loadcase/initial_stresses.h"
#include "loadcase/initial_velocity.h"
#include "loadcase/loadcase.h"
#include "loadcase/nodal_load_constant.h"
#include "loadcase/nodal_load_linear.h"
#include "loadcase/nodal_load_sin.h"
// #include "loadcase/uniaxial_load.h"
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
#include "material/spring_contact.h"
#include "material/spring_elastic.h"
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
// #include "reorderer/king.h"
// #include "reorderer/minimum_degree_ordering.h"
#include "reorderer/reorderer.h"
#include "reorderer/reverse_cuthill_mckee.h"
#include "reorderer/reverse_sloan.h"
#include "soe/band_linear_soe.h"
#include "soe/full_linear_soe.h"
#include "soe/symm_linear_soe.h"
#include "tracker/node_tracker.h"
#include "tracker/material_tracker.h"

/******************************************************************************
* Static variables
******************************************************************************/
static Domain*   pD = 0;
static Analysis* pA = 0;
static int currentLC = 0;
static int current_group = 0;

PyObject* buildList(const Vector& v) {
  static PyObject* pyList;
  pyList = PyList_New(v.size());
  for (int i = 0;i < v.size();i++)
    PyList_SetItem(pyList, i, PyFloat_FromDouble(v[i]));
  return pyList;
}

PyObject* buildList(const Matrix& m) {
  static PyObject* pyList;
  pyList = PyList_New(m.get_rows());
  for (int i = 0; i < m.get_rows(); i++) {
    PyObject* pyRow = PyList_New(m.get_cols());
    for (int j = 0;j < m.get_cols();j++)
      PyList_SetItem(pyRow, j, PyFloat_FromDouble(m(i, j)));
    PyList_SetItem(pyList, i, pyRow);
  }
  return pyList;
}

static PyObject* buildDict(std::istream& s) {
  static PyObject* pyDict;
  pyDict = PyDict_New();
  char name[128];
  int tag;
  s >> name;  // Read label, i.e. NODE
  s >> name;  // Read first tag name
  while (strcmp(name, "END")) {
    s >> tag;
    PyObject* pyKey = PyString_FromString(name);
    if (tag == 1000) {
      int n;
      s >> n;
      PyDict_SetItem(pyDict, pyKey, PyInt_FromLong(n));
    } else if (tag == 1020) {
      double d;
      s >> d;
      PyDict_SetItem(pyDict, pyKey, PyFloat_FromDouble(d));
    } else if (tag == 1100) {
      int n;
      s >> n;
      Vector v(n);
      for (int i = 0;i < n;i++) s >> v[i];
      PyDict_SetItem(pyDict, pyKey, buildList(v));
    } else if (tag == 1200) {
      int rows, cols;
      s >> rows;
      s >> cols;
      Matrix m(rows, cols, 0.);
      for (int i = 0;i < rows*cols;i++) s >> m.get_data()[i];
      PyDict_SetItem(pyDict, pyKey, buildList(m));
    } else {
      throw SException("[nemesis:%d] %s", 9999, "Internal error: Unknown tag.");
    }
    s >> name;
  }
  return pyDict;
}

PyObject* buildList(std::istream& s) {
  static PyObject* pyList;
  char name[128];
  int tag;
  s >> name;  // TRACKER
  s >> name;  // steps
  s >> tag;   // 1000
  int steps;
  s >> steps;
  pyList = PyList_New(steps);
  for (int i = 0; i < steps; i++) {
    double lambda, time;
    s >> name;  // lambda
    s >> tag;   // 1010
    s >> lambda;
    s >> name;  // time
    s >> tag;   // 1010
    s >> time;
    s >> name;  // data
    s >> tag;   // 1020
    PyObject* pyRow = PyList_New(3);
    PyList_SetItem(pyRow, 0, PyFloat_FromDouble(lambda));
    PyList_SetItem(pyRow, 1, PyFloat_FromDouble(time));
    PyList_SetItem(pyRow, 2, buildDict(s));
    PyList_SetItem(pyList, i, pyRow);
  }
  return pyList;
}
/******************************************************************************
* Database Commands
******************************************************************************/
static PyObject* pyDatabase_SQLite(PyObject* /*self*/, PyObject* args) {
  const char* s;
  if (!PyArg_ParseTuple(args, "s", &s))  return NULL;
  try {
    Database* pDB = new SQLiteDatabase(s);
    pD->set_database(pDB);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
    return Py_None;
}


static PyObject* pyDatabase_Store(PyObject* /*self*/, PyObject *args) {
  const char* s;
  if (!PyArg_ParseTuple(args, "s", &s))  return NULL;
  try {
    /// @todo Implement this method
    throw SException("[nemesis:%d] %s", 1110, "Unimplemented method");
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
   return Py_None;
}


static PyObject* pyDatabase_Close(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, ""))  return NULL;
  try {
    if (pD->get_database() == 0)
      throw SException("[nemesis:%d] %s", 1110, "No database set.");
    pD->get_database()->closeDB();
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyDatabase_ExportToVtk(PyObject* /*self*/, PyObject* args) {
  const char* s;
  if (!PyArg_ParseTuple(args, "s", &s))  return NULL;
  try {
    if (pD->get_database() == 0)
      throw SException("[nemesis:%d] %s", 1110, "No database set yet.");
    pD->get_database()->exportToVtk(s);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef DatabaseMethods[] =  {
  {"SQLite",  pyDatabase_SQLite,
    METH_VARARGS, "Create and utilize an SQLite database."},
  {"close", pyDatabase_Close,
    METH_VARARGS, "Close current database."},
  {"store", pyDatabase_Store,
    METH_VARARGS, "Store domain's state to the database."},
  {"exportToVtk", pyDatabase_ExportToVtk,
    METH_VARARGS, "Export a stored domain's state to vtk file format."},
  {NULL, NULL, 0, NULL}
};
/******************************************************************************
* Domain Commands
******************************************************************************/
static PyObject* pyDomain_Dim(PyObject* /*self*/, PyObject* args) {
  int n;
  if (!PyArg_ParseTuple(args, "i", &n))  return NULL;
  try {
    pD->set_dim(n);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyDomain_PlaneStress(PyObject* /*self*/, PyObject* args) {
  double t = 1.0;
  if (!PyArg_ParseTuple(args, "|d", &t)) return NULL;
  try {
    pD->set_tag(TAG_DOMAIN_PLANE_STRESS);
    pD->set_fac(t);
    pD->set_dim(2);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyDomain_PlaneStrain(PyObject* /*self*/, PyObject* args) {
  double t = 1.0;
  if (!PyArg_ParseTuple(args, "|d", &t)) return NULL;
  try {
    pD->set_tag(TAG_DOMAIN_PLANE_STRAIN);
    pD->set_fac(t);
    pD->set_dim(2);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyDomain_Axisymmetric(PyObject* /*self*/, PyObject* args) {
  double t = 1.0;
  if (!PyArg_ParseTuple(args, "|d", &t)) return NULL;
  try {
    pD->set_tag(TAG_DOMAIN_AXISYMMETRIC);
    pD->set_fac(t);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyDomain_RayleighDamping(PyObject* /*self*/, PyObject* args) {
  double aK = 0;
  double bM = 0;
  if (!PyArg_ParseTuple(args, "dd", &aK, &bM))  return NULL;
  Vector facs(2);
  facs[0]=aK;
  facs[1]=bM;
  pD->set_Rayleigh_factors(facs);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyDomain_Clear(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, ""))  return NULL;
  pD->clear();
  pA->clear();
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyDomain_State(PyObject* /*self*/, PyObject* args) {
  double facD;
  if (!PyArg_ParseTuple(args, "d", &facD)) return NULL;
  pD->state(facD);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyDomain_EigenValues(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  return buildList(pD->get_eigen_values());
}
static PyObject* pyDomain_Type(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  return Py_BuildValue("i", pD->get_tag());
}
static PyObject* pyDomain_Gravity(PyObject* /*self*/, PyObject* args) {
  double g, xG, yG, zG;
  if (!PyArg_ParseTuple(args, "d(ddd)", &g, &xG, &yG, &zG)) return NULL;
  try {
    pD->set_gravity(g, xG, yG, zG);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyDomain_GetState(PyObject* /*self*/, PyObject* args) {
  const char* s;
  if (!PyArg_ParseTuple(args, "")) return NULL;
  try {
    s = pD->get_state();
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  return Py_BuildValue("s", s);
}
static PyMethodDef DomainMethods[] =  {
  {"dim",   pyDomain_Dim,
    METH_VARARGS, "Set domain dimensions."},
  {"planeStress",   pyDomain_PlaneStress,
    METH_VARARGS, "Define a plain stress domain."},
  {"planeStrain",   pyDomain_PlaneStrain,
    METH_VARARGS, "Define a plain strain domain."},
  {"axisymmetric",    pyDomain_Axisymmetric,
    METH_VARARGS, "Define an axisymmetric domain."},
  {"clear", pyDomain_Clear,
    METH_VARARGS, "Clear the domain."},
  {"state", pyDomain_State,
    METH_VARARGS, "Set domain state."},
  {"RayleighDamping", pyDomain_RayleighDamping,
    METH_VARARGS, "Set Rayleigh damping factors."},
  {"eigenvalues",   pyDomain_EigenValues,
    METH_VARARGS, "Return the eigenvalues of the domain."},
  {"type",          pyDomain_Type,
    METH_VARARGS, "Domain type."},
  {"gravity",       pyDomain_Gravity,
    METH_VARARGS, "Set gravity direction."},
  {"get_state",     pyDomain_GetState,
    METH_VARARGS, "Return the domain's state."},
  {NULL, NULL, 0, NULL}
};
/******************************************************************************
* Node Commands
******************************************************************************/
static PyObject* pyNode_Add(PyObject* /*self*/, PyObject* args) {
  int id;
  double x1, x2 = 0, x3 = 0;
  if (!PyArg_ParseTuple(args, "id|dd", &id, &x1, &x2, &x3)) return NULL;
  try {
    Node* pNode = new Node(id, x1, x2, x3);
    pD->add(pD->get_nodes(), pNode);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyNode_Fix(PyObject* /*self*/, PyObject* args) {
  int node_id, dof;
  double c = 0.;
  if (!PyArg_ParseTuple(args, "ii|d", &node_id, &dof, &c))  return NULL;
  try {
    Node* node = pD->get<Node>(pD->get_nodes(), node_id);
    Constraint* constraint = new Constraint();
    constraint->set_cdof(node, dof, 1.0);
    constraint->set_val(c);
    pD->add(pD->get_constraints(), constraint);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyNode_Data(PyObject* /*self*/, PyObject* args) {
  int id;
  // define an output string stream
  if (!PyArg_ParseTuple(args, "i", &id)) return NULL;
  std::ostringstream os;
  try {
    pD->get<Node>(pD->get_nodes(), id)->Save(&os);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  static string tmp;
  tmp = os.str();
  return Py_BuildValue("s", tmp.c_str());
}


// static PyObject* pyNode_Path(PyObject* /*self*/, PyObject* args) {
//  int id;
//  const char* s;
//  if (!PyArg_ParseTuple(args, "i", &id)) return NULL;
//  try {
//    s = pD->get<Node>(pD->get_nodes(), id)->get_tracker()->data();
//  } catch(SException e) {
//    PyErr_SetString(PyExc_StandardError, e.what());
//    return NULL;
//  }
//  return Py_BuildValue("s", s);
// }

static PyMethodDef NodeMethods[] =  {
  {"add",   pyNode_Add,
    METH_VARARGS, "Adds a node to the domain."},
  {"fix",   pyNode_Fix,
    METH_VARARGS, "Fix a nodal degree of freedom."},
  {"data",  pyNode_Data,
    METH_VARARGS, "Returns a tuple containing nodal data."},
//  {"track", pyNode_Track,
//    METH_VARARGS, "Add a Tracker to the Node."},
//  {"path",  pyNode_Path,
//    METH_VARARGS, "Return the data recorded to the tracker."},
  {NULL, NULL, 0, NULL}
};

/*******************************************************************************
* CrossSection Commands
*******************************************************************************/
static PyObject* pySection_Rect(PyObject* /*self*/, PyObject* args) {
  int id;
  double w;
  double h;
  if (!PyArg_ParseTuple(args, "idd", &id, &w, &h)) return NULL;
  CrossSection* pSection = new RectangularCrossSection(id, w, h);
  pD->add(pD->get_sections(), pSection);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pySection_UserDefined(PyObject* /*self*/, PyObject* args) {
  int id;
  double A = 0, As2 = 0, As3 = 0, J1 = 0, J2 = 0, J3 = 0, h2 = 0, h3 = 0;
  if (!PyArg_ParseTuple(args, "id|ddddddd",
                               &id, &A, &As2, &As3,
                               &J1, &J2, &J3, &h2, &h3)) {
  return NULL;
  }
  CrossSection* pSection=
    new UserDefinedCrossSection(id, A, As2, As3, J1, J2, J3, h2, h3);
  pD->add(pD->get_sections(), pSection);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pySection_Rem(PyObject* /*self*/, PyObject* /*args*/) {
  return Py_None;
}
static PyObject* pySection_Info(PyObject* /*self*/, PyObject* args) {
  int id;
    if (!PyArg_ParseTuple(args, "i", &id)) return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef SectionMethods[] =  {
  {"user",  pySection_UserDefined,
    METH_VARARGS, "Adds a user defined section to the domain."},
  {"rect",  pySection_Rect,
    METH_VARARGS, "Adds a rectangular section to the domain."},
  {"rem",   pySection_Rem,
    METH_VARARGS, "Removes a section from the domain."},
  {"info",  pySection_Info,
    METH_VARARGS, "Prints section info."},
  {NULL, NULL, 0, NULL}
};

/*******************************************************************************
* Material Commands
*******************************************************************************/
static PyObject* pyMaterial_SDof(PyObject* /*self*/, PyObject* args) {
  int id;
  double E, rho = 0;
    if (!PyArg_ParseTuple(args, "id|d", &id, &E, &rho))return NULL;
  Material* pMaterial = new SDofMaterial(id, E, rho);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyMaterial_SpringElastic(PyObject* /*self*/, PyObject* args) {
  int id;
  double Kn, Ks2 = 0., Ks3 = 0.;
  if (!PyArg_ParseTuple(args, "id|dd", &id, &Kn, &Ks2, &Ks3)) {
    return NULL;
  }
  int dim = pD->get_dim();
  Material* pMaterial = new SpringElastic(id, dim, Kn, Ks2, Ks3);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_SpringContact(PyObject* /*self*/, PyObject* args) {
  int id;
  double Kn, Ks, mu, gap;
  if (!PyArg_ParseTuple(args, "idddd", &id, &Kn, &Ks, &mu, &gap)) {
    return NULL;
  }
  int dim = pD->get_dim();
  Material* pMaterial = new SpringContact(id, dim, Kn, Ks, mu, gap);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyMaterial_UniaxialElastic(PyObject* /*self*/,
                                            PyObject* args) {
  int id;
  double E = 0, nu = 0, rho = 0, aT = 0;
  if (!PyArg_ParseTuple(args, "id|ddd", &id, &E, &nu, &rho, &aT)) {
    return NULL;
  }
  Material* pMaterial = new UniaxialElastic(id, E, nu, rho, aT);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_UniaxialElastoPlastic(PyObject* /*self*/,
                                                  PyObject* args) {
  int id;
  double E = 0, nu = 0, rho = 0, aT = 0, sy = 0, Hiso = 0, Hkin = 0, eta = 0;
  if (!PyArg_ParseTuple(args, "iddddddd|d",
                        &id, &E, &nu, &rho, &aT, &sy, &Hiso, &Hkin, &eta)) {
    return NULL;
  }
  Material* pMaterial = new UniaxialElastoPlastic(id, E, nu, rho, aT,
                                                  sy, Hiso, Hkin, eta);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_UniaxialGap(PyObject* /*self*/, PyObject* args) {
  int id;
  double E = 0, nu = 0, rho = 0, aT = 0, sy = 0, gap = 0;
  if (!PyArg_ParseTuple(args, "idddddd", &id, &E, &nu, &rho, &aT, &sy, &gap))
    return NULL;
  Material* pMaterial = new UniaxialGap(id, E, nu, rho, aT, sy, gap);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_Elastic(PyObject* /*self*/, PyObject *args) {
  int id;
  double E = 0, nu = 0, rho = 0, aT = 0;
  double kx = 0., ky = 0., kz = 0.;
  if (!PyArg_ParseTuple(args, "idd|ddddd",
                        &id, &E, &nu, &rho, &aT, &kx, &ky, &kz)) {
    return NULL;
  }
  Material* pMaterial;
  if (pD->get_tag() == TAG_DOMAIN_PLANE_STRESS)
    pMaterial = new PlaneStress(id, E, nu, rho, aT);
  else
    pMaterial = new MultiaxialElastic(id, E, nu, rho, aT, kx, ky, kz);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_DuncanChang(PyObject* /*self*/, PyObject* args) {
  int id;
  double E, nu, c, phi, m, Rf, pa, rho = 0, aT = 0;
  if (!PyArg_ParseTuple(args, "iddddddd|dd",
                        &id, &E, &nu, &c, &phi, &m, &Rf, &pa, &rho, &aT)) {
    return NULL;
  }
  Material* pMaterial = new DuncanChang(id, E, nu, c, phi, m, Rf, pa, rho, aT);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_UniaxialCyclic(PyObject* /*self*/, PyObject* args) {
  int id;
  double E, nu, rho, aT, tmax, gmax;
  if (!PyArg_ParseTuple(args, "idddddd", &id, &E, &nu, &rho, &aT, &tmax, &gmax))
    return NULL;
  Material* pMaterial = new UniaxialCyclic(id, E, nu, rho, aT, tmax, gmax);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_VonMises(PyObject* /*self*/, PyObject* args) {
  int id, elasticId;
  double s0, K;
  if (!PyArg_ParseTuple(args, "iidd", &id, &elasticId, &s0, &K))
    return NULL;
  MultiaxialMaterial* elastic;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), elasticId);
    elastic = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Material* material = new VonMises(id, elastic, s0, K);
  pD->add(pD->get_materials(), material);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_MohrCoulomb(PyObject* /*self*/, PyObject* args) {
  int id, elasticId;
  double c, phi, alpha;
  if (!PyArg_ParseTuple(args, "iiddd", &id, &elasticId, &c, &phi, &alpha))
    return NULL;
  MultiaxialMaterial* elastic;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), elasticId);
    elastic = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Material* pMaterial = new MohrCoulomb(id, elastic, c, phi, alpha);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_Tresca(PyObject* /*self*/, PyObject* args) {
  int id, elasticId;
  double cu;
  double kx = 0., ky = 0., kz = 0.;
  if (!PyArg_ParseTuple(args, "iid|ddd", &id, &elasticId, &cu, &kx, &ky, &kz))
    return NULL;
  MultiaxialMaterial* elastic;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), elasticId);
    elastic = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Material* pMaterial = new Tresca(id, elastic, cu, kx, ky, kz);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_HoekBrown(PyObject* /*self*/, PyObject* args) {
  int id, elasticId;
  double si, sp, mb, mbb, alpha;
  if (!PyArg_ParseTuple(args, "iiddddd",
                        &id, &elasticId, &si, &sp, &mb, &mbb, &alpha)) {
    return NULL;
  }
  MultiaxialMaterial* elastic;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), elasticId);
    elastic = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Material* pMaterial = new HoekBrown(id, elastic, si, sp, mb, mbb, alpha);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_Creep(PyObject* /*self*/, PyObject* args) {
  int id, elasticId;
  double A, n, k;
  if (!PyArg_ParseTuple(args, "iiddd", &id, &elasticId, &A, &n, &k))
    return NULL;
  MultiaxialMaterial* elastic;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), elasticId);
    elastic = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Material* pMaterial = new Creep(id, elastic, A, n, k);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_DruckerPrager(PyObject* /*self*/, PyObject* args) {
  int id, elasticId, type;
  double c, phi, psi, T;
  if (!PyArg_ParseTuple(args, "iiidddd",
                        &id, &elasticId, &type, &c, &phi, &psi, &T)) {
    return NULL;
  }
  MultiaxialMaterial* elastic;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), elasticId);
    elastic = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Material* pMaterial = new DruckerPrager(id, elastic, type, c, phi, psi, T);
  pD->add(pD->get_materials(), pMaterial);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_DruckerPragerNew(PyObject* /*self*/,
                                             PyObject* args) {
  int id, elasticId;
  double c, phi, psi, Kphi, Kc, T;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iidddddd",
                        &id, &elasticId, &c, &phi, &psi, &Kc, &Kphi, &T)) {
    return NULL;
  }
  MultiaxialMaterial* elastic;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), elasticId);
    elastic = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create material and add it to domain
  Material* pMaterial = new DruckerPragerNew(id, elastic,
                                             c, phi, psi, Kc, Kphi, T);
  pD->add(pD->get_materials(), pMaterial);
  // Increase reference count and return
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_DruckerPragerNew2(PyObject* /*self*/,
                                              PyObject* args) {
  int id, elasticId;
  double c, phi, psi, Kphi, Kc, T;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iidddddd",
                        &id, &elasticId, &c, &phi, &psi, &Kc, &Kphi, &T)) {
    return NULL;
  }
  MultiaxialMaterial* elastic;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), elasticId);
    elastic = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create material and add it to domain
  Material* pMaterial = new DruckerPragerNew2(id, elastic,
                                              c, phi, psi, Kc, Kphi, T);
  pD->add(pD->get_materials(), pMaterial);
  // Increase reference count and return
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_DruckerPragerNew3(PyObject* /*self*/,
                                               PyObject* args) {
  int id, elasticId;
  double c, phi, psi, Kphi, Kc, T;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iidddddd",
                        &id, &elasticId, &c, &phi, &psi, &Kc, &Kphi, &T)) {
    return NULL;
  }
  MultiaxialMaterial* elastic;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), elasticId);
    elastic = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create material and add it to domain
  Material* mat = new DruckerPragerNew3(id, elastic,
                                        c, phi, psi, Kc, Kphi, T);
  pD->add(pD->get_materials(), mat);
  // Increase reference count and return
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_ModifiedCamClay(PyObject* /*self*/,
                                            PyObject* args) {
  int id, elasticId;
  double M, po, kappa, lambda;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iidddd",
                        &id, &elasticId, &M, &po, &kappa, &lambda))
    return NULL;
  MultiaxialMaterial* elastic;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), elasticId);
    elastic = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create material and add it to domain
  Material* mat = new ModifiedCamClay(id, elastic, M, po, kappa, lambda);
  pD->add(pD->get_materials(), mat);
  // Increase reference count and return
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyMaterial_LadeDuncan(PyObject* /*self*/, PyObject* args) {
  int id, elasticId;
  double K;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iid", &id, &elasticId, &K)) {
    return NULL;
  }
  MultiaxialMaterial* elastic;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), elasticId);
    elastic = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create material and add it to domain
  Material* pMaterial = new LadeDuncan(id, elastic, K);
  pD->add(pD->get_materials(), pMaterial);
  // Increase reference count and return
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef MaterialMethods[] =  {
  {"sdof",            pyMaterial_SDof,
    METH_VARARGS, "Define a single dof material."},
  {"springElastic",       pyMaterial_SpringElastic,
    METH_VARARGS, "Define a spring elastic material."},
  {"springContact",       pyMaterial_SpringContact,
    METH_VARARGS, "Define a spring contact material."},
  {"uniElastic",          pyMaterial_UniaxialElastic,
    METH_VARARGS, "Define a uniaxial elastic material."},
  {"uniElastoPlastic",      pyMaterial_UniaxialElastoPlastic,
    METH_VARARGS, "Define a uniaxial elastoplastic material."},
  {"uniCyclic",         pyMaterial_UniaxialCyclic,
    METH_VARARGS, "Define a cyclic material based on a backbone curve."},
  {"uniGap",            pyMaterial_UniaxialGap,
    METH_VARARGS, "Define a uniaxial elastoplastic material with a gap."},
  {"elastic",           pyMaterial_Elastic,
    METH_VARARGS, "Define a multiaxial elastic material."},
  {"DuncanChang",         pyMaterial_DuncanChang,
    METH_VARARGS, "Define a Duncan-Chang material."},
  {"vonMises",          pyMaterial_VonMises,
    METH_VARARGS, "Define a von-Mises type material."},
  {"MohrCoulomb",         pyMaterial_MohrCoulomb,
    METH_VARARGS, "Define a Mohr-Coulomb type material."},
  {"Tresca",            pyMaterial_Tresca,
    METH_VARARGS, "Define a Tresca type material."},
  {"HoekBrown",           pyMaterial_HoekBrown,
    METH_VARARGS, "Define a Hoek Brown type material."},
  {"DruckerPrager",       pyMaterial_DruckerPrager,
    METH_VARARGS, "Define a Drucker-Prager type material."},
  {"DruckerPragerNew",      pyMaterial_DruckerPragerNew,
    METH_VARARGS, "Define a Drucker-Prager type material."},
  {"DruckerPragerNew2",     pyMaterial_DruckerPragerNew2,
    METH_VARARGS, "Define a Drucker-Prager type material."},
  {"DruckerPragerNew3",     pyMaterial_DruckerPragerNew3,
    METH_VARARGS, "Define a Drucker-Prager type material."},
  {"modifiedCamClay",       pyMaterial_ModifiedCamClay,
    METH_VARARGS, "Define a Modified Cam-Clay type material."},
  {"LadeDuncan",          pyMaterial_LadeDuncan,
    METH_VARARGS, "Define a Lade-Duncan type material."},
  {"creep",           pyMaterial_Creep,
    METH_VARARGS, "Define a creep material."},
  {NULL, NULL, 0, NULL}
};

/*******************************************************************************
* Element commands
*******************************************************************************/
/**
 *
 */
static PyObject* pyElement_Spring(PyObject* /*self*/, PyObject* args) {
  int id, na, nb, mat_id, dim;
  double xp1 = 1., xp2 = 0., xp3 = 0., yp1 = 0., yp2 = 1., yp3 = 0.;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iiiii|(ddd)(ddd)", &id, &na, &nb, &mat_id, &dim,
    &xp1, &xp2, &xp3, &yp1, &yp2, &yp3)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(2);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), na);
    nodes[1] = pD->get<Node>(pD->get_nodes(), nb);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  SpringMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat_id);
    /// @todo Check material before casting
    material = static_cast<SpringMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element and add it to domain
  try {
    Spring* spring = new Spring(id, nodes, material, dim,
                                xp1, xp2, xp3, yp1, yp2, yp3);
    pD->add(pD->get_elements(), spring);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(spring);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 *
 */
static PyObject* pyElement_Bar2s(PyObject* /*self*/, PyObject* args) {
  int id, n1, n2, mat, s1, s2=-9999;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iiiii|i", &id, &n1, &n2, &mat, &s1, &s2)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(2);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  UniaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<UniaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get section pointers from domain
  CrossSection* sec1;
  CrossSection* sec2;
  try {
    sec1 = pD->get<CrossSection>(pD->get_sections(), s1);
    if (s2 == -9999) {
      sec2 = sec1;
    } else {
      sec2 = pD->get<CrossSection>(pD->get_sections(), s2);
    }
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get domain dimension
  int dim = pD->get_dim();
  // Create element and add it to domain
  try {
    Bar2s* elem = new Bar2s(id, nodes, material, sec1, sec2, dim);
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 *
 */
static PyObject* pyElement_Bar2t(PyObject* /*self*/, PyObject* args) {
  int id, n1, n2, mat, s1, s2=-9999;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iiiii|i", &id, &n1, &n2, &mat, &s1, &s2)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(2);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  UniaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<UniaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get section pointers from domain
  CrossSection* sec1;
  CrossSection* sec2;
  try {
    sec1 = pD->get<CrossSection>(pD->get_sections(), s1);
    if (s2 == -9999) {
      sec2 = sec1;
    } else {
      sec2 = pD->get<CrossSection>(pD->get_sections(), s2);
    }
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get domain dimension
  int dim = pD->get_dim();
  // Create element and add it to domain
  try {
    Bar2t* elem = new Bar2t(id, nodes, material, sec1, sec2, dim);
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Create a 2-node Euler Beam.
 */
static PyObject* pyElement_Beam2e(PyObject* /*self*/, PyObject* args) {
  int id, n1, n2, mat, sec;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iiiii", &id, &n1, &n2, &mat, &sec)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(2);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  UniaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<UniaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get section pointers from domain
  CrossSection* section;
  try {
    section = pD->get<CrossSection>(pD->get_sections(), sec);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element and add it to domain
  try {
    Beam2e* elem = new Beam2e(id, nodes, material, section);
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Create a 2-node Timoshenko Beam.
 */
static PyObject* pyElement_Beam2t(PyObject* /*self*/, PyObject* args) {
  int id, n1, n2, mat, sec, rule = 1;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iiiii|i", &id, &n1, &n2, &mat, &sec, &rule)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(2);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  UniaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<UniaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get section pointers from domain
  CrossSection* section;
  try {
    section = pD->get<CrossSection>(pD->get_sections(), sec);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element and add it to domain
  try {
    Timoshenko2d* elem = new Timoshenko2d(id, nodes, material, section, rule);
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Create a 3-node Timoshenko Beam.
 */
static PyObject* pyElement_Beam3t(PyObject* /*self*/, PyObject* args) {
  int id, n1, n2, n3, mat, sec, rule = 2;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iiiiii|i", &id, &n1, &n2, &n3,
                                          &mat, &sec, &rule)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(2);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  UniaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<UniaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get section pointers from domain
  CrossSection* section;
  try {
    section = pD->get<CrossSection>(pD->get_sections(), sec);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element and add it to domain
  try {
    Timoshenko2d* elem = new Timoshenko2d(id, nodes, material, section, rule);
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Create an 8-node brick.
 *
 */
static PyObject* pyElement_Brick8d(PyObject* /*self*/, PyObject* args) {
  // Local variables
  int id, n1, n2, n3, n4, n5, n6, n7, n8, mat;
  // Get data
  if (!PyArg_ParseTuple(args, "iiiiiiiiii",
                        &id, &n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8, &mat)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(8);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
    nodes[3] = pD->get<Node>(pD->get_nodes(), n4);
    nodes[4] = pD->get<Node>(pD->get_nodes(), n5);
    nodes[5] = pD->get<Node>(pD->get_nodes(), n6);
    nodes[6] = pD->get<Node>(pD->get_nodes(), n7);
    nodes[7] = pD->get<Node>(pD->get_nodes(), n8);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  MultiaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element
  Brick8d* elem = new Brick8d(id, nodes, material);
  // Add element to the domain/group
  try {
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Create an 8-node brick.
 *
 */
static PyObject* pyElement_Brick8b(PyObject* /*self*/, PyObject* args) {
  // Local variables
  int id, n1, n2, n3, n4, n5, n6, n7, n8, mat;
  // Get data
  if (!PyArg_ParseTuple(args, "iiiiiiiiii",
                        &id, &n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8, &mat)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(8);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
    nodes[3] = pD->get<Node>(pD->get_nodes(), n4);
    nodes[4] = pD->get<Node>(pD->get_nodes(), n5);
    nodes[5] = pD->get<Node>(pD->get_nodes(), n6);
    nodes[6] = pD->get<Node>(pD->get_nodes(), n7);
    nodes[7] = pD->get<Node>(pD->get_nodes(), n8);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  MultiaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element
  Brick8b* elem = new Brick8b(id, nodes, material);
  // Add element to the domain/group
  try {
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Create an 8-node brick.
 *
 */
static PyObject* pyElement_Brick8i(PyObject* /*self*/, PyObject* args) {
  // Local variables
  int id, n1, n2, n3, n4, n5, n6, n7, n8, mat;
  // Get data
  if (!PyArg_ParseTuple(args, "iiiiiiiiii",
                        &id, &n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8, &mat)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(8);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
    nodes[3] = pD->get<Node>(pD->get_nodes(), n4);
    nodes[4] = pD->get<Node>(pD->get_nodes(), n5);
    nodes[5] = pD->get<Node>(pD->get_nodes(), n6);
    nodes[6] = pD->get<Node>(pD->get_nodes(), n7);
    nodes[7] = pD->get<Node>(pD->get_nodes(), n8);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  MultiaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element
  Brick8i* elem = new Brick8i(id, nodes, material);
  // Add element to the domain/group
  try {
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Create a single dof element.
 */
static PyObject* pyElement_SDOF(PyObject* /*self*/, PyObject* args) {
  int id, node, dof, mat;
  // Check for consistent input
  if (!PyArg_ParseTuple(args, "iiii", &id, &node, &dof, &mat)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(2);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), node);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  SDofMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<SDofMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element and add it to domain
  try {
    SDofElement* elem = new SDofElement(id, nodes, dof, material);
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Add a new quad (plane stress/plane strain/axisymmetric to the domain)
 *
 */
static PyObject* pyElement_Quad4d(PyObject* /*self*/, PyObject* args) {
  // Local variables
  int id, n1, n2, n3, n4, mat_id;
  double thickness = 1.0;
  // Get data
  if (!PyArg_ParseTuple(args, "iiiiii|d",
                        &id, &n1, &n2, &n3, &n4, &mat_id, &thickness)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(4);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
    nodes[3] = pD->get<Node>(pD->get_nodes(), n4);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  MultiaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat_id);
    /// @todo Check material before casting
    material = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element
  Quad4* quad4;
  try {
    switch (pD->get_tag()) {
      case TAG_DOMAIN_PLANE_STRAIN:
      case TAG_DOMAIN_PLANE_STRESS:
        quad4 = new Quad4DispPlain(id, nodes, material, thickness);
        break;
      case TAG_DOMAIN_AXISYMMETRIC:
        quad4 = new Quad4DispAxisymmetric(id, nodes, material, thickness);
        break;
      default:
        throw SException("[nemesis:%d] %s", 1110,
          "A quad4d can be used only in plane strain/stress/axisymmetry.");
    }
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Add element to the domain/group
  try {
    pD->add(pD->get_elements(), quad4);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(quad4);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Add a new quad (plane stress/plane strain/axisymmetric) to the domain
 *
 */
static PyObject* pyElement_Quad4Test(PyObject* /*self*/, PyObject* args) {
  // Local variables
  int id, n1, n2, n3, n4, mat_id;
  double thickness = 1.0;
  // Get data
  if (!PyArg_ParseTuple(args, "iiiiii|d",
                        &id, &n1, &n2, &n3, &n4, &mat_id, &thickness)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(4);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
    nodes[3] = pD->get<Node>(pD->get_nodes(), n4);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  MultiaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat_id);
    /// @todo Check material before casting
    material = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element
  Quad4* quad4;
  try {
    switch (pD->get_tag()) {
      case TAG_DOMAIN_PLANE_STRAIN:
      case TAG_DOMAIN_PLANE_STRESS:
        quad4 = new Quad4d(id, nodes, material, thickness, false);
        break;
      case TAG_DOMAIN_AXISYMMETRIC:
        quad4 = new Quad4d(id, nodes, material, thickness, true);
        break;
      default:
        throw SException("[nemesis:%d] %s", 1110,
          "A quad4d can be used only in plane strain/stress/axisymmetry.");
    }
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Add element to the domain/group
  try {
    pD->add(pD->get_elements(), quad4);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(quad4);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Add a new b-bar quad (plane stress/plane strain/axisymmetric to the domain)
 *
 */
static PyObject* pyElement_Quad4b(PyObject* /*self*/, PyObject* args) {
  // Local variables
  int id, n1, n2, n3, n4, mat_id;
  double thickness = 1.0;
  // Get data
  if (!PyArg_ParseTuple(args, "iiiiii|d",
                        &id, &n1, &n2, &n3, &n4, &mat_id, &thickness)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(4);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
    nodes[3] = pD->get<Node>(pD->get_nodes(), n4);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  MultiaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat_id);
    /// @todo Check material before casting
    material = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element
  Quad4* quad4;
  try {
    switch (pD->get_tag()) {
      case TAG_DOMAIN_PLANE_STRAIN:
      case TAG_DOMAIN_PLANE_STRESS:
        quad4 = new Quad4b(id, nodes, material, thickness, false);
        break;
      case TAG_DOMAIN_AXISYMMETRIC:
        quad4 = new Quad4b(id, nodes, material, thickness, true);
        break;
      default:
        throw SException("[nemesis:%d] %s", 1110,
          "A quad4d can be used only in plane strain/stress/axisymmetry.");
    }
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Add element to the domain/group
  try {
    pD->add(pD->get_elements(), quad4);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(quad4);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Add a new quad (incompatible nodes) to the domain
 *
 */
static PyObject* pyElement_Quad4i(PyObject* /*self*/, PyObject* args) {
  // Local variables
  int id, n1, n2, n3, n4, mat_id;
  double thickness = 1.0;
  // Get data
  if (!PyArg_ParseTuple(args, "iiiiii|d",
                        &id, &n1, &n2, &n3, &n4, &mat_id, &thickness)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(4);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
    nodes[3] = pD->get<Node>(pD->get_nodes(), n4);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  MultiaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat_id);
    /// @todo Check material before casting
    material = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element
  Quad4* quad4;
  try {
    switch (pD->get_tag()) {
      /// @todo:  check
      case TAG_DOMAIN_PLANE_STRAIN:
      case TAG_DOMAIN_PLANE_STRESS:
        quad4 = new Quad4i(id, nodes, material, thickness);
        break;
      case TAG_DOMAIN_AXISYMMETRIC:
      default:
        throw SException("[nemesis:%d] %s", 1110,
          "A quad4d can be used only in plane strain/stress.");
    }
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Add element to the domain/group
  try {
    pD->add(pD->get_elements(), quad4);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(quad4);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Add a new quad (incompatible nodes) to the domain
 *
 */
static PyObject* pyElement_Quad4e(PyObject* /*self*/, PyObject* args) {
  // Local variables
  int id, n1, n2, n3, n4, mat_id;
  double thickness = 1.0;
  // Get data
  if (!PyArg_ParseTuple(args, "iiiiii|d",
                        &id, &n1, &n2, &n3, &n4, &mat_id, &thickness)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(4);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
    nodes[3] = pD->get<Node>(pD->get_nodes(), n4);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  MultiaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat_id);
    /// @todo Check material before casting
    material = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element
  Quad4* quad4;
  try {
    switch (pD->get_tag()) {
      /// @todo:  check
      case TAG_DOMAIN_PLANE_STRAIN:
      case TAG_DOMAIN_PLANE_STRESS:
        quad4 = new Quad4e(id, nodes, material, thickness);
        break;
      case TAG_DOMAIN_AXISYMMETRIC:
      default:
        throw SException("[nemesis:%d] %s", 1110,
          "A quad4e can be used only in plane strain/stress.");
    }
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Add element to the domain/group
  try {
    pD->add(pD->get_elements(), quad4);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(quad4);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Create a new 3-node triangle
 * @todo Check for the axisymmetric case 
 */
static PyObject* pyElement_Triangle3d(PyObject* /*self*/, PyObject* args) {
  // Local variables
  int id, n1, n2, n3, mat;
  double thickness = 1.0;
  // Get data
  if (!PyArg_ParseTuple(args, "iiiii|d",
                        &id, &n1, &n2, &n3, &mat, &thickness)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(3);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  MultiaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element
  Triangle3* elem = new Triangle3(id, nodes, material, thickness);
  // Add element to the domain/group
  try {
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Create a new 6-node triangle
 * @todo Check for the axisymmetric case 
 */
static PyObject* pyElement_Triangle6d(PyObject* /*self*/, PyObject* args) {
  // Local variables
  int id, n1, n2, n3, n4, n5, n6, mat;
  double thickness = 1.0;
  // Get data
  if (!PyArg_ParseTuple(args, "iiiiiiii|d",
                        &id, &n1, &n2, &n3, &n4, &n5, &n6, &mat, &thickness)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(3);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
    nodes[3] = pD->get<Node>(pD->get_nodes(), n4);
    nodes[4] = pD->get<Node>(pD->get_nodes(), n5);
    nodes[5] = pD->get<Node>(pD->get_nodes(), n6);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  MultiaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element
  Triangle6* elem = new Triangle6(id, nodes, material, thickness);
  // Add element to the domain/group
  try {
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Create an 8-node brick.
 *
 */
static PyObject* pyElement_Tetrahedron4d(PyObject* /*self*/, PyObject* args) {
  // Local variables
  int id, n1, n2, n3, n4, mat;
  // Get data
  if (!PyArg_ParseTuple(args, "iiiiii", &id, &n1, &n2, &n3, &n4, &mat)) {
    return NULL;
  }
  // Get node pointers from domain
  std::vector<Node*> nodes = std::vector<Node*>(8);
  try {
    nodes[0] = pD->get<Node>(pD->get_nodes(), n1);
    nodes[1] = pD->get<Node>(pD->get_nodes(), n2);
    nodes[2] = pD->get<Node>(pD->get_nodes(), n3);
    nodes[3] = pD->get<Node>(pD->get_nodes(), n4);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Get material pointer from domain
  MultiaxialMaterial* material;
  try {
    Material* m = pD->get<Material>(pD->get_materials(), mat);
    /// @todo Check material before casting
    material = static_cast<MultiaxialMaterial*>(m);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  // Create element
  Tetrahedron4Disp* elem = new Tetrahedron4Disp(id, nodes, material);
  // Add element to the domain/group
  try {
    pD->add(pD->get_elements(), elem);
    Group* group = pD->get<Group>(pD->get_groups(), current_group);
    group->AddElement(elem);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyElement_Data(PyObject* /*self*/, PyObject* args) {
  int id;
  // define an output string stream
  if (!PyArg_ParseTuple(args, "i", &id)) return NULL;
  std::ostringstream os;
  try {
    pD->get<Element>(pD->get_elements(), id)->Save(&os);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  static string tmp;
  tmp = os.str();
  return Py_BuildValue("s", tmp.c_str());
}
static PyMethodDef ElementMethods[] = {
  {"spring",      pyElement_Spring,
    METH_VARARGS, "Define a spring element."},
  {"bar",       pyElement_Bar2s,
    METH_VARARGS, "An alias for bar2s."},
  {"bar2s",     pyElement_Bar2s,
    METH_VARARGS, "Defines a simple 1d/2d/3d geometrically linear bar."},
  {"bar2t",     pyElement_Bar2t,
    METH_VARARGS, "Defines a 1d/2d/3d bar (Total Lagrangian formulation)."},
  {"beam2e",      pyElement_Beam2e,
    METH_VARARGS, "Defines a 2-Node Euler-Bernulli beam."},
  {"beam2t",      pyElement_Beam2t,
    METH_VARARGS, "Defines a 2-Node Timoshenko beam."},
  {"beam3t",      pyElement_Beam3t,
    METH_VARARGS, "Defines a 3-Node Timoshenko beam."},
  {"sdof",      pyElement_SDOF,
    METH_VARARGS, "Define a single dof element."},
  {"quad4d",      pyElement_Quad4d,
    METH_VARARGS, "Define a 4-Noded standard displacement quad."},
  {"quad4t",      pyElement_Quad4Test,
    METH_VARARGS, "Define a 4-Noded standard displacement quad."},
  {"quad4b",      pyElement_Quad4b,
    METH_VARARGS, "Define a 4-Noded B-Bar quad."},
  {"quad4i",      pyElement_Quad4i,
    METH_VARARGS, "Define a 4-Noded non-conforming quad."},
  {"quad4e",      pyElement_Quad4e,
    METH_VARARGS, "Define a 4-Noded enhanced assumed strain quad."},
  {"tria3d",        pyElement_Triangle3d,
    METH_VARARGS, "Define a 3-Noded constant strain triangle."},
  {"tria6d",        pyElement_Triangle6d,
    METH_VARARGS, "Define a 6-Noded linear strain triangle."},
  {"tetra4d",     pyElement_Tetrahedron4d,
    METH_VARARGS, "Define a 4-Noded constant strain tetrahedron."},
  {"brick8d", pyElement_Brick8d,
    METH_VARARGS, "Define a 8-Noded standard displacement brick."},
  {"brick8b", pyElement_Brick8b,
    METH_VARARGS, "Define a 8-Noded Bbar brick."},
  {"brick8i", pyElement_Brick8i,
    METH_VARARGS, "Define a 8-Noded non-conforming brick."},
  {"data",      pyElement_Data,
    METH_VARARGS, "Access to the element data."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* Constraint commands
******************************************************************************/
static PyObject* pyConstraint_Set(PyObject* /*self*/, PyObject* args) {
  int node_id, dof;
  double a, c;
  if (!PyArg_ParseTuple(args, "d(iid)", &c, &node_id, &dof, &a))
    return NULL;
  try {
    Node* node = pD->get<Node>(pD->get_nodes(), node_id);
    Constraint* constraint = new Constraint;
    constraint->set_cdof(node, dof, a);
    constraint->set_val(c);
    pD->add(pD->get_constraints(), constraint);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
/// @todo make this one more general
static PyObject* pyConstraint_twoDofs(PyObject* /*self*/, PyObject* args) {
  int node_id1, dof1, node_id2, dof2;
  double a1, a2, c;
  if (!PyArg_ParseTuple(args, "(iid)(iid)d", &node_id1, &dof1, &a1,
                                             &node_id2, &dof2, &a2, &c))
    return NULL;
  try {
    Node* node1 = pD->get<Node>(pD->get_nodes(), node_id1);
    Node* node2 = pD->get<Node>(pD->get_nodes(), node_id2);
    Constraint* constraint = new Constraint;
    constraint->set_cdof(node1, dof1, a1);
    constraint->set_cdof(node2, dof2, a2);
    constraint->set_val(c);
    pD->add(pD->get_constraints(), constraint);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyConstraint_Linear(PyObject* /*self*/, PyObject* args) {
  double c0;
  PyObject* cList;
  if (!PyArg_ParseTuple(args, "dO", &c0, &cList)) return NULL;
  Constraint* constraint = new Constraint;
  constraint->set_val(c0);

  // Check if list ok
  if (!PyList_Check(cList)) {
    PyErr_SetString(PyExc_StandardError, "Expected a list.");
    return NULL;
  }
  // Start reading list
  for (int i = 0; i < PyList_Size(cList); i++) {
    // Get tuple
    PyObject* cTuple = PyList_GetItem(cList, i);
    // Check if tuple ok
    if (!PyTuple_Check(cTuple) || (PyTuple_Size(cTuple) !=3)) {
      PyErr_SetString(PyExc_StandardError, "Expected a tuple of size 3.");
      return NULL;
    }
    // Get arguments and check if arguments ok
    PyObject* tmpNode = PyTuple_GetItem(cTuple, 0);
    PyObject* tmpDof =PyTuple_GetItem(cTuple, 1);
    PyObject* tmpCi  =PyTuple_GetItem(cTuple, 2);
    if (!PyInt_Check(tmpNode)) {
      PyErr_SetString(PyExc_StandardError, "Expected (INT, int, float).");
      return NULL;
    }
    if (!PyInt_Check(tmpDof)) {
      PyErr_SetString(PyExc_StandardError, "Expected (int, INT, float).");
      return NULL;
    }
    if (!PyFloat_Check(tmpCi)) {
      PyErr_SetString(PyExc_StandardError, "Expected (int, int, FLOAT).");
      return NULL;
    }
    // Get arguments
    int node_id = PyInt_AsLong(tmpNode);
    int dof     = PyInt_AsLong(tmpDof);
    double ci   = PyFloat_AsDouble(tmpCi);
    // Build triple
    try {
      Node* node = pD->get<Node>(pD->get_nodes(), node_id);
      constraint->set_cdof(node, dof, ci);
    } catch(SException e) {
      PyErr_SetString(PyExc_StandardError, e.what());
      return NULL;
    }
  }
  // Add to domain
  pD->add(pD->get_constraints(), constraint);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef ConstraintMethods[] =  {
  {"set",             pyConstraint_Set,
    METH_VARARGS, "Define a single/multiple dof constraint."},
  {"twoDofs",           pyConstraint_twoDofs,
    METH_VARARGS, "Define a two dof constraint."},
  {"linear",            pyConstraint_Linear,
    METH_VARARGS, "Define a linear constraint."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* Group commands
******************************************************************************/
static PyObject* pyGroup_Define(PyObject* /*self*/, PyObject* args) {
  int id;
  if (!PyArg_ParseTuple(args, "i", &id)) return NULL;
  try {
    Group* group = new Group(id);
    pD->add(pD->get_groups(), group);
    current_group = id;
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyGroup_State(PyObject* /*self*/, PyObject* args) {
  if (currentLC <= 0) {
    PyErr_SetString(PyExc_StandardError, "No LoadCase yet defined.");
    return NULL;
  }
  int id;
  int active;
  double facK = 1.;
  double facS = 1.;
  double facG = 1.;
  double facP = 1.;
  if (!PyArg_ParseTuple(args, "ii|ddddd", &id, &active,
                        &facK, &facS, &facG, &facP))
    return NULL;
  try {
    Group* group = pD->get<Group>(pD->get_groups(), id);
    GroupState* gs = new GroupState(group, active, facK, facS, facG, facP);
    pD->get<LoadCase>(pD->get_loadcases(), currentLC)->AddGroupState(gs);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef GroupMethods[] =  {
  {"define",  pyGroup_Define,
    METH_VARARGS, "Define and set current group."},
  {"state", pyGroup_State,
    METH_VARARGS, "Define group state for current laodcase."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* Sensitivity commands
******************************************************************************/
static PyObject* pySens_Elem(PyObject* /*self*/, PyObject* args) {
  int id, parameter;
  if (!PyArg_ParseTuple(args, "ii", &id, &parameter)) {
    return NULL;
  }
  try {
    Element* e = pD->get<Element>(pD->get_elements(), id);
    LoadCase* lc = pD->get<LoadCase>(pD->get_loadcases(), currentLC);
    lc->AddSensitivityParameter(new ElementSensitivityParameter(e, parameter));
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 *
 */
static PyMethodDef SensitivityMethods[] =  {
  {"elem",  pySens_Elem,
    METH_VARARGS, "Define element sensitivity parameter."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* Load commands
******************************************************************************/
static PyObject* pyLoad_NodeConstant(PyObject* /*self*/, PyObject* args) {
  if (currentLC <= 0) {
    PyErr_SetString(PyExc_StandardError, "No LoadCase yet defined.");
    return NULL;
  }
  int node_id, dof;
  double value;
  if (!PyArg_ParseTuple(args, "iid", &node_id, &dof, &value))
    return NULL;
  try {
    Node* node = pD->get<Node>(pD->get_nodes(), node_id);
    Load* load = new NodalLoadConstant(node, dof, value);
    pD->get<LoadCase>(pD->get_loadcases(), currentLC)->AddLoad(load);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyLoad_Linear(PyObject* /*self*/, PyObject* args) {
  if (currentLC <= 0) {
    PyErr_SetString(PyExc_StandardError, "No LoadCase yet defined.");
    return NULL;
  }
  int node_id, dof;
  double initial_value, grad;
  if (!PyArg_ParseTuple(args, "iidd", &node_id, &dof, &initial_value, &grad))
    return NULL;
  try {
    Node* node = pD->get<Node>(pD->get_nodes(), node_id);
    Load* load = new NodalLoadLinear(node, dof, initial_value, grad);
    pD->get<LoadCase>(pD->get_loadcases(), currentLC)->AddLoad(load);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyLoad_Sin(PyObject* /*self*/, PyObject* args) {
  if (currentLC <= 0) {
    PyErr_SetString(PyExc_StandardError, "No LoadCase yet defined.");
    return NULL;
  }
  int node_id, dof;
  double a, omega, phi;
  if (!PyArg_ParseTuple(args, "iiddd", &node_id, &dof, &a, &omega, &phi))
    return NULL;
  try {
    Node* node = pD->get<Node>(pD->get_nodes(), node_id);
    Load* load = new NodalLoadSin(node, dof, a, omega, phi);
    pD->get<LoadCase>(pD->get_loadcases(), currentLC)->AddLoad(load);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef LoadMethods[] =  {
  {"node",    pyLoad_NodeConstant,
    METH_VARARGS, "Define a constant nodal load."},
  {"linear",    pyLoad_Linear,
    METH_VARARGS, "Define a linear in time nodal load."},
  {"sin",     pyLoad_Sin,
    METH_VARARGS, "Define a sinus nodal load."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* Ground Motion commands
******************************************************************************/
static PyObject* pyGroundMotion_File(PyObject* /*self*/, PyObject* args) {
  if (currentLC <= 0)  /// @todo remove exception to lc
    throw SException("[nemesis:%d] %s", 1110, "No LoadCase yet defined.");
  int dof;
  const char* filename;
  double dt;
  double scale = 1.0;
  if (!PyArg_ParseTuple(args, "isd|d", &dof, &filename, &dt, &scale))
    return NULL;
  try {
    std::ifstream data(filename);
    const std::map<int, Element*>* elements = &(pD->get_elements());
    Load* load = new GroundMotionFile(elements, dof, data, dt, scale);
    data.close();
    pD->get<LoadCase>(pD->get_loadcases(), currentLC)->AddLoad(load);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyGroundMotion_Sin(PyObject* /*self*/, PyObject* args) {
  if (currentLC <= 0)  /// @todo remove exception to lc
    throw SException("[nemesis:%d] %s", 1110, "No LoadCase yet defined.");
  int dof;
  double a, omega, phi = 0.;
  if (!PyArg_ParseTuple(args, "idd|d", &dof, &a, &omega, &phi))
    return NULL;
  try {
    const std::map<int, Element*>* elements = &(pD->get_elements());
    Load* load = new GroundMotionSin(elements, dof, a, omega, phi);
    pD->get<LoadCase>(pD->get_loadcases(), currentLC)->AddLoad(load);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef GroundMotionMethods[] =  {
  {"file",  pyGroundMotion_File,
    METH_VARARGS, "Define a uniform excitation through a file."},
  {"sin",   pyGroundMotion_Sin,
    METH_VARARGS, "Define a sinus uniform excitation."},
  {NULL, NULL, 0, NULL}
};

/*******************************************************************************
* Initial Conditions commands
*******************************************************************************/
static PyObject* pyInitialConditions_Disp(PyObject* /*self*/, PyObject* args) {
  if (currentLC <= 0) {
    PyErr_SetString(PyExc_StandardError, "No LoadCase yet defined.");
    return NULL;
  }
  int node_id, dof;
  double disp;
  if (!PyArg_ParseTuple(args, "iid", &node_id, &dof, &disp))
    return NULL;
  try {
    Node* node = pD->get<Node>(pD->get_nodes(), node_id);
    InitialCondition* ic = new InitialDisplacement(node, dof, disp);
    pD->get<LoadCase>(pD->get_loadcases(), currentLC)->AddInitialCondition(ic);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyInitialConditions_Velc(PyObject* /*self*/, PyObject* args) {
  if (currentLC <= 0) {
    PyErr_SetString(PyExc_StandardError, "No LoadCase yet defined.");
    return NULL;
  }
  int node_id, dof;
  double velocity;
  if (!PyArg_ParseTuple(args, "iid", &node_id, &dof, &velocity))
    return NULL;
  try {
    Node* node = pD->get<Node>(pD->get_nodes(), node_id);
    InitialCondition* ic = new InitialVelocity(node, dof, velocity);
    pD->get<LoadCase>(pD->get_loadcases(), currentLC)->AddInitialCondition(ic);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyInitialConditions_Stresses(PyObject* /*self*/,
                                              PyObject* args) {
  if (currentLC <= 0) {
    PyErr_SetString(PyExc_StandardError, "No LoadCase yet defined.");
    return NULL;
  }
  int id, dir;
  double h1, s1, h2, s2, K0;
  if (!PyArg_ParseTuple(args, "iiddddd", &id, &dir, &h1, &s1, &h2, &s2, &K0))
    return NULL;
  try {
    Group* group = pD->get<Group>(pD->get_groups(), id);
    InitialStresses* is = new InitialStresses(group, dir, h1, s1, h2, s2, K0);
    pD->get<LoadCase>(pD->get_loadcases(), currentLC)->AddInitialCondition(is);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef InitialConditionsMethods[] =  {
  {"disp",    pyInitialConditions_Disp,
    METH_VARARGS, "Set initial displacement to nodal dof."},
  {"velc",    pyInitialConditions_Velc,
    METH_VARARGS, "Set initial velocity to nodal dof."},
  {"stresses",    pyInitialConditions_Stresses,
    METH_VARARGS, "Set initial stresses to a group."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* LoadCase commands
******************************************************************************/
static PyObject* pyLC_Define(PyObject* /*self*/, PyObject* args) {
  int id;
  const char* label = 0;
  if (!PyArg_ParseTuple(args, "i|s", &id, &label)) {
    return NULL;
  }
  LoadCase* loadcase;
  /// @todo: check this better
  if (label != 0) {
    loadcase = new LoadCase(id, label);
  } else {
    loadcase = new LoadCase(id, "default");
  }
  pD->add(pD->get_loadcases(), loadcase);
  currentLC = id;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyLC_Info(PyObject* /*self*/, PyObject* args) {
  int id;
  if (!PyArg_ParseTuple(args, "i", &id)) {
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef LCMethods[] =  {
  {"define",  pyLC_Define,
    METH_VARARGS, "Define a loadcase."},
  {"info",  pyLC_Info,
    METH_VARARGS, "Print loadcase info."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* Analysis commands
******************************************************************************/
static PyObject* pyAnalysis_Static(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  AnalysisType* pType = new StaticAnalysis();
  pA->set_analysis_type(pType);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyAnalysis_Transient(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  AnalysisType* pType = new TransientAnalysis();
  pA->set_analysis_type(pType);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyAnalysis_Eigen(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  AnalysisType* pType = new EigenAnalysis();
  pA->set_analysis_type(pType);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyAnalysis_Sensitivity(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  AnalysisType* pType = new SensitivityStaticAnalysis();
  try {
    pA->set_analysis_type(pType);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyAnalysis_Run(PyObject* /*self*/, PyObject* args) {
  int id, steps, ret;
  if (!PyArg_ParseTuple(args, "ii", &id, &steps)) return NULL;
  try {
    LoadCase* loadcase = pD->get<LoadCase>(pD->get_loadcases(), id);
    ret = pA->analyze(loadcase, steps);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef AnalysisMethods[] =  {
  {"static",  pyAnalysis_Static,
    METH_VARARGS, "Define a static analysis."},
  {"transient", pyAnalysis_Transient,
    METH_VARARGS, "Define a transient analysis."},
  {"eigen",   pyAnalysis_Eigen,
    METH_VARARGS, "Define a eigeinvalue analysis."},
  {"sensitivity", pyAnalysis_Sensitivity,
    METH_VARARGS, "Define a sensitivity analysis."},
  {"run", pyAnalysis_Run,
    METH_VARARGS, "Run given analysis."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* Imposer commands
******************************************************************************/
static PyObject* pyImposer_Eliminatation(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  Imposer* pImposer = new EliminationImposer();
  pA->set_imposer(pImposer);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyImposer_Lagrange(PyObject* /*self*/, PyObject* args) {
  double a;
  if (!PyArg_ParseTuple(args, "", &a)) return NULL;
  Imposer* pImposer = new LagrangeImposer();
  pA->set_imposer(pImposer);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyImposer_Penalty(PyObject* /*self*/, PyObject* args) {
  double a;
  if (!PyArg_ParseTuple(args, "d", &a)) return NULL;
  Imposer* pImposer = new PenaltyImposer(a);
  pA->set_imposer(pImposer);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef ImposerMethods[] =  {
  {"elimination", pyImposer_Eliminatation,
    METH_VARARGS, "Use the elimination method to impose constraints."},
  {"lagrange",  pyImposer_Lagrange,
    METH_VARARGS, "Use the Langange multipliers method to impose constraints."},
  {"penalty", pyImposer_Penalty,
    METH_VARARGS, "Use the penalty method to impose constraints."},
  {NULL, NULL, 0, NULL}
};

/*******************************************************************************
* Control commands
*******************************************************************************/
static PyObject* pyControl_Load(PyObject* /*self*/, PyObject* args) {
  double DL0, minDL = 1., maxDL = 1., n = 0.5, dt = 0.;
  int Id = 1;
  if (!PyArg_ParseTuple(args, "d|ddidd", &DL0, &minDL, &maxDL, &Id, &n, &dt))
    return NULL;
  Control* pControl = new LoadControl(DL0, minDL, maxDL, Id, n, dt);
  pA->set_control(pControl);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyControl_ArcLength(PyObject* /*self*/, PyObject* args) {
  double DL0, minDL = 1., maxDL = 1., n = 0.5, dt = 0.;
  int Id = 1;
  if (!PyArg_ParseTuple(args, "d|ddidd", &DL0, &minDL, &maxDL, &Id, &n, &dt))
    return NULL;
  Control* pControl = new ArcLengthSpherical(DL0, minDL, maxDL, Id, n, dt);
  pA->set_control(pControl);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyControl_ArcLengthUNP(PyObject* /*self*/, PyObject* args) {
  double DL0, minDL = 1., maxDL = 1., n = 0.5, dt = 0.;
  int Id = 1;
  if (!PyArg_ParseTuple(args, "d|ddidd", &DL0, &minDL, &maxDL, &Id, &n, &dt))
    return NULL;
  Control* pControl = new ArcLengthUNP(DL0, minDL, maxDL, Id, n, dt);
  pA->set_control(pControl);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyControl_Disp(PyObject* /*self*/, PyObject* args) {
  int node, dof;
  double Du0, minDu = 1., maxDu = 1., n = 0.5, dt = 0.;
  int Id = 1;
  if (!PyArg_ParseTuple(args, "iid|ddidd",
                        &node, &dof, &Du0, &minDu, &maxDu, &Id, &n, &dt)) {
    return NULL;
  }
  Control* c = new DisplacementControl(node, dof,
                                       Du0, minDu, maxDu, Id, n, dt);
  pA->set_control(c);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyControl_Newmark(PyObject* /*self*/, PyObject* args) {
  double beta, gamma, dt;
  if (!PyArg_ParseTuple(args, "ddd", &beta, &gamma, &dt)) return NULL;
  Control* pControl = new Newmark(beta, gamma, dt);
  pA->set_control(pControl);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef ControlMethods[] = {
  {"load",      pyControl_Load,
    METH_VARARGS, "Use load control for the analysis."},
  {"arcLength",   pyControl_ArcLength,
    METH_VARARGS, "Use arc-length control for the analysis."},
  {"arcLengthUNP",    pyControl_ArcLengthUNP,
    METH_VARARGS, "Use arc-length control for the analysis."},
  {"disp",      pyControl_Disp,
    METH_VARARGS, "Use arc-length control for the analysis."},
  {"Newmark",     pyControl_Newmark,
    METH_VARARGS, "Use Newmark control for the analysis."},
  {NULL, NULL, 0, NULL}
};

/*******************************************************************************
* Algorithm commands
*******************************************************************************/
static PyObject* pyAlgorithm_Linear(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  Algorithm* pAlgorithm = new LinearAlgorithm();
  pA->set_algorithm(pAlgorithm);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyAlgorithm_BFGS(PyObject* /*self*/, PyObject* args) {
  int m;
  int lsearch = 0;
  double etaMin = 0.1, etaMax = 1.0, rTol = 0.8;
  int maxIter = 10;
  if (!PyArg_ParseTuple(args, "i|idddi",
                        &m, &lsearch, &etaMin, &etaMax, &rTol, &maxIter)) {
    return NULL;
  }
  Algorithm* pAlgorithm;
  if (lsearch == 0)  pAlgorithm = new BFGS(m);
  else      pAlgorithm = new BFGS(m, etaMin, etaMax, rTol, maxIter);
  pA->set_algorithm(pAlgorithm);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyAlgorithm_FullNewtonRaphson(PyObject* /*self*/,
                                               PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  Algorithm* pAlgorithm = new NewtonRaphsonFull();
  pA->set_algorithm(pAlgorithm);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyAlgorithm_InitialNewtonRaphson(PyObject* /*self*/,
                                                  PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  Algorithm* pAlgorithm = new NewtonRaphsonInitial();
  pA->set_algorithm(pAlgorithm);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
pyAlgorithm_ModifiedNewtonRaphson(PyObject* /*self*/,
                                  PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  Algorithm* pAlgorithm = new NewtonRaphsonModified();
  pA->set_algorithm(pAlgorithm);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef AlgorithmMethods[] =  {
  {"linear",      pyAlgorithm_Linear,
    METH_VARARGS, "Use linear algorithm."},
  {"BFGS",      pyAlgorithm_BFGS,
    METH_VARARGS, "Use BFGS rank-two update algorithm."},
  {"fNR",       pyAlgorithm_FullNewtonRaphson,
    METH_VARARGS, "Use full Newton-Raphson algorithm."},
  {"mNR",       pyAlgorithm_ModifiedNewtonRaphson,
    METH_VARARGS, "Use modified Newton-Raphson algorithm."},
  {"iNR",       pyAlgorithm_InitialNewtonRaphson,
    METH_VARARGS, "Use Newton-Raphson algorithm with initial stiffness."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* Convergence commands
******************************************************************************/
static PyObject* pyConvergence_Set(PyObject* /*self*/, PyObject* args) {
  int maxIter;
  double tolRabs;
  double tolRrel = 1.0;
  double tolWrel = 1.0;
  if (!PyArg_ParseTuple(args, "id|dd",
    &maxIter, &tolRabs, &tolRrel, &tolWrel)) return NULL;
  pA->get_convergence_norm()->set_check(maxIter, tolRabs, tolRrel, tolWrel);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef ConvergenceMethods[] =  {
  {"set", pyConvergence_Set,
    METH_VARARGS, "Set convergence check's properties."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* SOE commands
******************************************************************************/
static PyObject* pySOE_FullLinearSOE(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  SOE* pSOE = new FullLinearSOE();
  pA->set_soe(pSOE);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pySOE_SymmLinearSOE(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  SOE* pSOE = new SymmLinearSOE();
  pA->set_soe(pSOE);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pySOE_BandLinearSOE(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  SOE* pSOE = new BandLinearSOE();
  pA->set_soe(pSOE);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef SOEMethods[] =  {
  {"full",  pySOE_FullLinearSOE,
    METH_VARARGS, "Use full storage scheme for the system of equations."},
  {"symm",  pySOE_SymmLinearSOE,
    METH_VARARGS, "Use symmetric storage scheme for the system of equations."},
  {"band",  pySOE_BandLinearSOE,
    METH_VARARGS, "Use band storage scheme for the system of equations."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* Reorder commands
******************************************************************************/
static PyObject* pyReorder_RCM(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  Reorderer* pReorderer = new ReverseCuthillMckee();
  pA->set_reorderer(pReorderer);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyReorder_FCM(PyObject* /*self*/, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  Reorderer* pReorderer = new ForwardCuthillMckee();
  pA->set_reorderer(pReorderer);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyReorder_FSloan(PyObject* /*self*/, PyObject* args) {
  double w1 = 0.5, w2 = 0.5;
  if (!PyArg_ParseTuple(args, "|dd", &w1, &w2)) return NULL;
  Reorderer* pReorderer = new ForwardSloan(w1, w2);
  pA->set_reorderer(pReorderer);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* pyReorder_RSloan(PyObject* /*self*/, PyObject* args) {
  double w1 = 0.5, w2 = 0.5;
  if (!PyArg_ParseTuple(args, "|dd", &w1, &w2)) return NULL;
  Reorderer* pReorderer = new ReverseSloan(w1, w2);
  pA->set_reorderer(pReorderer);
  Py_INCREF(Py_None);
  return Py_None;
}
/*
static PyObject* pyReorder_King(PyObject* self, PyObject* args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  Reorderer* pReorderer = new King();
  pA->set_reorderer(pReorderer);
  Py_INCREF(Py_None);
  return Py_None;
}
*/
static PyMethodDef ReorderMethods[] =  {
  {"rCM", pyReorder_RCM,
    METH_VARARGS, "Use reverse Cuthill-Mckee to reorder nodal numbering."},
  {"fCM", pyReorder_FCM,
    METH_VARARGS, "Use forward Cuthill-Mckee to reorder nodal numbering."},
  {"rSloan",  pyReorder_RSloan,
    METH_VARARGS, "Use reverse Sloan to reorder nodal numbering."},
  {"fSloan",  pyReorder_FSloan,
    METH_VARARGS, "Use forward Sloan to reorder nodal numbering."},
//  {"mdo", pyReorder_MinimumDegreeOrdering,
//    METH_VARARGS, "Use Minimum Degree Ordering to reorder nodal numbering."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* Tracker commands
******************************************************************************/
static PyObject* pyTracker_Node(PyObject* /*self*/, PyObject* args) {
  // Parse data
  int id, node_id;
  if (!PyArg_ParseTuple(args, "ii", &id, &node_id)) {
    return NULL;
  }
  // Get node pointer from the domain and then
  // create tracker and add it to the domain
  try {
    Node* node = pD->get<Node>(pD->get_nodes(), node_id);
    Tracker* tracker = new NodeTracker(id, node);
    pD->add<Tracker>(pD->get_trackers(), tracker);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyTracker_MatPoint(PyObject* /*self*/, PyObject* args) {
  // Parse data
  int id, elem_id, mat_id;
  if (!PyArg_ParseTuple(args, "iii", &id, &elem_id, &mat_id)) {
    return NULL;
  }
  /// @todo Add material tracker
  // Get material pointer from and then
  // create tracker and add it to the domain
  // try {
    // Element* elem = pD->get<Element>(pD->get_elements(), id);
    // Tracker* tracker = new MaterialTracker(id, material);
    // pD->add<Tracker>(pD->get_trackers(), tracker);
  // } catch(SException e) {
  //  PyErr_SetString(PyExc_StandardError, e.what());
  //  return NULL;
  // }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* pyTracker_Data(PyObject* /*self*/, PyObject* args) {
  int id;
  // define an output string stream
  if (!PyArg_ParseTuple(args, "i", &id)) return NULL;
  std::ostringstream os;
  try {
    pD->get<Tracker>(pD->get_trackers(), id)->Save(&os);
  } catch(SException e) {
    PyErr_SetString(PyExc_StandardError, e.what());
    return NULL;
  }
  static string tmp;
  tmp = os.str();
  return Py_BuildValue("s", tmp.c_str());
}

static PyMethodDef TrackerMethods[] =  {
  {"node",      pyTracker_Node,
    METH_VARARGS, "Add a new tracker to given node."},
  {"matpoint",      pyTracker_MatPoint,
    METH_VARARGS, "Add a new tracker to given material point."},
  {"data",      pyTracker_Data,
    METH_VARARGS, "Get tracker data for a given step."},
  {NULL, NULL, 0, NULL}
};

/******************************************************************************
* Log commands
******************************************************************************/
PyObject* log_CaptureStdout(PyObject* /*self*/, PyObject* pArgs) {
  char* LogStr = NULL;
  if (!PyArg_ParseTuple(pArgs, "s", &LogStr)) return NULL;
  printf("%s", LogStr);
  // Simply using printf to do the real work.
  // You could also write it to a .log file or whatever...
  // MessageBox(NULL, LogStr...
  // WriteFile(hFile, LogStr...
  Py_INCREF(Py_None);
  return Py_None;
}

// Notice we have STDERR too.
PyObject* log_CaptureStderr(PyObject* /*self*/, PyObject* pArgs) {
  char* LogStr = NULL;
  if (!PyArg_ParseTuple(pArgs, "s", &LogStr)) return NULL;
  printf("%s", LogStr);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef logMethods[] = {
  {"CaptureStdout", log_CaptureStdout, METH_VARARGS, "Logs stdout"},
  {"CaptureStderr", log_CaptureStderr, METH_VARARGS, "Logs stderr"},
  {NULL, NULL, 0, NULL}
};
/******************************************************************************
* PyParser Methods
******************************************************************************/
/**
 * Default constructor.
 */
PyParser::PyParser() {
  pD=&D;
  pA=&A;
}
/**
 * Destructor.
 */
PyParser::~PyParser() {
}
/**
* Parse interactively.
*/
int PyParser::parse() {
  Py_Initialize();
  this->initModules();
  FILE* fp = stdin;
  PyRun_InteractiveLoop(fp, "<Solver-Interactive>");
  Py_Finalize();
  return 0;
}
/**
* Parse file.
*/
int PyParser::parse(char* filename) {
  // Check if the extension is correct
  int len = strlen(filename);
  if ((  filename[len-4] != '.'
    ||filename[len-3] != 's'
    ||filename[len-2] != 'l'
    ||filename[len-1] != 'v')) {
    printf("solver: Files passed with '-p' should have an '.slv' extension.");
    return 0;
  }
  // Check if file exists
  std::fstream tmpFile;
  tmpFile.open(filename, std::ios::binary|std::ios::in);
  if (tmpFile.fail()) {
    tmpFile.clear();
    tmpFile.close();
    printf("solver: File not found!");
    return 0;
  }
  tmpFile.close();
  // Now it is ok to parse the file
  Py_Initialize();
  this->initModules();
  char mode[] = "r";  // This removes deprecated warning.
  PyObject* PyFileObject = PyFile_FromString(filename, mode);
  PyRun_SimpleFile(PyFile_AsFile(PyFileObject), filename);
  Py_DECREF(PyFileObject);
  Py_Finalize();
  return 0;
}
int PyParser::initModules() {
  Py_InitModule("db", DatabaseMethods);
  Py_InitModule("domain", DomainMethods);
  Py_InitModule("node", NodeMethods);
  Py_InitModule("section", SectionMethods);
  Py_InitModule("material", MaterialMethods);
  Py_InitModule("element", ElementMethods);
  Py_InitModule("constraint", ConstraintMethods);
  Py_InitModule("load", LoadMethods);
  Py_InitModule("groundMotion", GroundMotionMethods);
  Py_InitModule("initial", InitialConditionsMethods);
  Py_InitModule("lc", LCMethods);
  Py_InitModule("analysis", AnalysisMethods);
  Py_InitModule("imposer", ImposerMethods);
  Py_InitModule("control", ControlMethods);
  Py_InitModule("algorithm", AlgorithmMethods);
  Py_InitModule("convergence", ConvergenceMethods);
  Py_InitModule("soe", SOEMethods);
  Py_InitModule("reorder", ReorderMethods);
  Py_InitModule("group", GroupMethods);
  Py_InitModule("sens", SensitivityMethods);
  Py_InitModule("tracker", TrackerMethods);
  PyRun_SimpleString("import db");
  PyRun_SimpleString("import domain");
  PyRun_SimpleString("import node");
  PyRun_SimpleString("import section");
  PyRun_SimpleString("import material");
  PyRun_SimpleString("import material");
  PyRun_SimpleString("import element");
  PyRun_SimpleString("import constraint");
  PyRun_SimpleString("import load");
  PyRun_SimpleString("import groundMotion");
  PyRun_SimpleString("import initial");
  PyRun_SimpleString("import lc");
  PyRun_SimpleString("import analysis");
  PyRun_SimpleString("import imposer");
  PyRun_SimpleString("import control");
  PyRun_SimpleString("import algorithm");
  PyRun_SimpleString("import convergence");
  PyRun_SimpleString("import soe");
  PyRun_SimpleString("import reorder");
  PyRun_SimpleString("import group");
  PyRun_SimpleString("import sens");
  PyRun_SimpleString("import tracker");

  Py_InitModule("log", logMethods);
  PyRun_SimpleString("import log\n"
                     "import sys\n"
                     "class StdoutCatcher:\n"
                     "\tdef write(self, str):\n"
                     "\t\tlog.CaptureStdout(str)\n"
                     "class StderrCatcher:\n"
                     "\tdef write(self, str):\n"
                     "\t\tlog.CaptureStderr(str)\n"
                     "sys.stdout = StdoutCatcher()\n"
                     "sys.stderr = StderrCatcher()\n");
  return 0;
}
