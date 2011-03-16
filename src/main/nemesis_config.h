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

#ifndef SRC_MAIN_NEMESIS_CONFIG_H_
#define SRC_MAIN_NEMESIS_CONFIG_H_

/*****************************************************************************
* nemesis Info
*****************************************************************************/
#define NEMESIS_VERSION  "0.9.07"

/*****************************************************************************
* Platform Info
*****************************************************************************/
#if defined(WIN32)
  #define NEMESIS_PLATFORM "win32"
#elif defined(LINUX)
  #define NEMESIS_PLATFORM "linux"
#else
  #define NEMESIS_PLATFORM "unknown"
#endif

/*****************************************************************************
* Compiler Info
*****************************************************************************/
#if defined(__GNUC__)
  #define NEMESIS_COMPILER_ID 1
  #define NEMESIS_COMPILER_NAME "gcc"
  #define NEMESIS_COMPILER_VERSION __VERSION__
#elif defined(_MSC_VER)
  #define NEMESIS_COMPILER_ID 2
  #define NEMESIS_COMPILER_NAME "MSC"
  #define NEMESIS_COMPILER_VERSION _MSC_VER
#else
  #define NEMESIS_COMPILER_ID 0
  #define NEMESIS_COMPILER_NAME "<Unknown>"
  #define NEMESIS_COMPILER_VERSION 0.0.0
#endif

/*****************************************************************************
* snsprintf
*****************************************************************************/
// #if _MSC_VER >= 1400  // VC 8.0 and later deprecate snprintf and _snprintf.
//   #define snprintf _snprintf_s
// #elif _MSC_VER
//   #define snprintf _snprintf
// #endif

#endif  // SRC_MAIN_NEMESIS_CONFIG_H_
