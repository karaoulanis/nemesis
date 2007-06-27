/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 2, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License along   *
*   with this program; if not, write to the Free Software Foundation, Inc.,   *
*   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.               *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//*****************************************************************************

enum FEObjectTag{
	TAG_NONE									= 0000,
	TAG_NOT_YET_SPECIFIED						= 0001,
	
	// Node Tags								
	TAG_NODE									= 1000,
	
	// CrossSection Tags						
	TAG_CROSSSECTION_USER_DEFINED				= 1301,
	TAG_CROSSSECTION_RECTANGLE					= 1302,
	
	// Material Tags							
	TAG_MATERIAL_UNIAXIAL_ELASTIC				= 1501,
	TAG_MATERIAL_UNIAXIAL_ELASTOPLASTIC			= 1502,
	TAG_MATERIAL_UNIAXIAL_CYCLIC				= 1503,
	TAG_MATERIAL_SINGLE_DOF						= 1551,
	TAG_MATERIAL_MULTIAXIAL_ELASTIC				= 1601,
	TAG_MATERIAL_MOHR_COULOMB					= 1603,
	TAG_MATERIAL_NON_MISES						= 1604,
	TAG_MATERIAL_DRUCKER_PRAGER_IN				= 1605,
	TAG_MATERIAL_DRUCKER_PRAGER_OUT				= 1606,
	TAG_MATERIAL_TRESCA							= 1607,

	// Element Tags
	TAG_ELEM_BAR_2D_GEOMETRICALLY_LINEAR		= 2001,
	TAG_ELEM_BAR_2D_TOTAL_LAGRANGIAN			= 2002,
	TAG_ELEM_BAR_3D_GEOMETRICALLY_LINEAR		= 2003,
	TAG_ELEM_BAR_3D_TOTAL_LAGRANGIAN			= 2004,
	TAG_ELEM_BEAM_2D_EULER						= 2005,
	TAG_ELEM_BEAM_3D_EULER						= 2006,
	TAG_ELEM_QUAD_4_STD_DISP					= 2007,
	TAG_ELEM_TRIANGLE_3_PRESSURE				= 2008,
	TAG_ELEM_TRIANGLE_6							= 2009,
	TAG_ELEM_BRICK_8_DISP						= 2010,
	TAG_ELEM_TETRHEDRON_4_DISP					= 2011,
	TAG_ELEM_BEAM_2D_TIMOSHENKO_2				= 2012,
	TAG_ELEM_BEAM_2D_TIMOSHENKO_3				= 2013,

	// MaterialItem Tags						
	TAG_MATERIAL_POINT							= 4000,

	// Loads									
	TAG_NODAL_LOAD_CONSTANT						= 8001,
	TAG_NODAL_LOAD_SIN   						= 8003,
	TAG_LOAD_BAR_FORCE_CONCENTRATED				= 8011,
	TAG_LOAD_BAR_FORCE_UNIFORM					= 8012,
	TAG_LOAD_BAR_FORCE_TRAPEZOIDAL				= 8013,
	TAG_LOAD_BEAM_FORCE_CONCENTRATED			= 8014,
	TAG_LOAD_BEAM_FORCE_UNIFORM					= 8015,
	TAG_LOAD_BEAM_FORCE_TRAPEZOIDAL				= 8016,
	TAG_LOAD_BEAM_MOMENT_CONCENTRATED			= 8017,
	TAG_LOAD_BEAM_MOMENT_UNIFORM				= 8018,
	TAG_LOAD_BEAM_MOMENT_TRAPEZOIDAL			= 8019,
	TAG_LOAD_TEMPERATURE_UNIFORM				= 8020,
	TAG_LOAD_TEMPERATURE_GRADIENT				= 8021,
	TAG_LOAD_GRAVITY							= 8022,

	// Initial Conditions
	TAG_INITIAL_DISPLACEMENT					= 8050,
	TAG_INITIAL_VELOCITY						= 8051,
	TAG_INITIAL_STRESSES						= 8052,

	// Analysis Tags							
	TAG_ANALYSIS_STATIC							= 9001,
	TAG_ANALYSIS_TRANSIENT						= 9002,
	TAG_IMPOSER_ELIMINATION						= 9051,
	TAG_IMPOSER_PENALTY							= 9052,
	TAG_IMPOSER_LAGRANGE						= 9053,
	TAG_CONTROL_LOAD							= 9101,
	TAG_CONTROL_DISPLACEMENT					= 9102,
	TAG_CONTROL_ARC_LENGTH_SPHERICAL			= 9103,
	TAG_CONTROL_ARC_LENGTH_UNP					= 9104,
	TAG_CONTROL_GENERALIZED_A					= 9151,
	TAG_CONTROL_HAMILTON_VELOCITIES				= 9152,
	TAG_CONTROL_HAMILTON_DISPLACEMENTS			= 9153,
	TAG_ALGORITHM_LINEAR						= 9201,
	TAG_ALGORITHM_NEWTON_RAPHSON_FULL			= 9202,
	TAG_ALGORITHM_NEWTON_RAPHSON_MODIFED		= 9203,
	TAG_ALGORITHM_NEWTON_RAPHSON_INITIAL		= 9204,
	TAG_ALGORITHM_NEWTON_RAPHSON_PERIODIC		= 9205,
	TAG_SOE_FULL_GENERIC_POSITIVE_DEFINE		= 9301,
	TAG_SOE_LINEAR_FULL							= 9302,
	TAG_SOE_LINEAR_SYMM							= 9303,
	TAG_SOE_LINEAR_BAND							= 9304,
	TAG_CONVERGENCE_NORM_EUCLIDEAN				= 9501,
	TAG_CONVERGENCE_NORM_MAXIMUM				= 9502,
	TAG_TRACKER									= 9601,
	TAG_REORDERER_NONE							= 9701,
	TAG_REORDERER_FORWARD_CUTHILL_MCKEE			= 9702,
	TAG_REORDERER_REVERSE_CUTHILL_MCKEE			= 9703,
	TAG_REORDERER_FORWARD_SLOAN					= 9704,
	TAG_REORDERER_REVERSE_SLOAN					= 9705,
//	TAG_REORDERER_MINIMUM_DEGREE_ORDERING		= 9706,
//	TAG_REORDERER_KING							= 9707,
	TAG_PACKET									= 9901
};
