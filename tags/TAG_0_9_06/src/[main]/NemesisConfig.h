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

/*****************************************************************************
* nemesis Info
*****************************************************************************/
#define _NEMESIS_VERSION	"0.9.06"

/*****************************************************************************
* Platform Info
*****************************************************************************/
#if defined(WIN32)
	#define _NEMESIS_PLATFORM "win32"
#elif defined(LINUX)
	#define _NEMESIS_PLATFORM "linux"
#else
	#define _NEMESIS_PLATFORM "unknown"
#endif

/*****************************************************************************
* Compiler Info
*****************************************************************************/
#if defined(__GNUC__)
	#define _NEMESIS_COMPILER_ID 1
	#define _NEMESIS_COMPILER_NAME "gcc"
	#define _NEMESIS_COMPILER_VERSION __VERSION__
#elif defined(_MSC_VER)
	#define _NEMESIS_COMPILER_ID 2
	#define _NEMESIS_COMPILER_NAME "MSC"
	#define _NEMESIS_COMPILER_VERSION _MSC_VER
#else
	#define _NEMESIS_COMPILER_ID 0
	#define _NEMESIS_COMPILER_NAME "<Unknown>"
	#define _NEMESIS_COMPILER_VERSION 0.0.0
#endif

/*****************************************************************************
* Memory leak detectors
*****************************************************************************/
#ifdef _MSC_VER
//	#include <vld.h>
//	#include <vldapi.h>
#endif
