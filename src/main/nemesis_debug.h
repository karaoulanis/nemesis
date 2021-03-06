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

#ifndef SRC_MAIN_NEMESIS_DEBUG_H_
#define SRC_MAIN_NEMESIS_DEBUG_H_

#include "numeric/matrix.h"
#include "numeric/vector.h"

void report(const double  d, const char* name="Noname",
            int total = 8, int decimal = 4);
void report(const Matrix& m, const char* name="Noname",
            int total = 8, int decimal = 4);
void report(const Vector& v,
            const char* name="Noname", int total = 8, int decimal = 4);

void add(Matrix* K, int row, int col, const Matrix& B1,
         const Matrix& C, const Matrix B2, double c1, double c0 = 0.);
void add(Vector* R, int row, const Matrix& BT,
         const Vector& V, double c1, double c0 = 0.);
void add2(Vector* R, int row, const Matrix& BT, const Vector& V,
          double c1, double c0 = 0.);

void add_BTCB(Matrix* K, int row, int col, const int* perm,
              const Matrix& B1, const Matrix& C, const Matrix B2,
              double c1, double c0 = 0.);
void add_BTv(Vector* R, int row, const int* perm,
             const Matrix& B, const Vector& v,
             double c1, double c0 = 0.);
void add_Bv(Vector* R, int row, const int* perm,
            const Matrix& B, const Vector& v,
            double c1, double c0 = 0.);

void spectralDecomposition(const Vector& s,
                           Vector& sP, Matrix& sV);
#endif  // SRC_MAIN_NEMESIS_DEBUG_H_
