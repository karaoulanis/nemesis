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

// C++ system files
#include <fstream>
#include <iostream>
#include <string>

// Project files (alphabetically)
#include "numeric/matrix.h"
#include "numeric/vector.h"

class LogFile {
 private:
  std::ofstream log;
 public:
  explicit LogFile(const char* name)       {log.open(name, std::ios_base::out);}
  ~LogFile()                                                      {log.close();}
  LogFile& operator<<(char c)                          {log << c; return *this;}
  LogFile& operator<<(unsigned char c)                 {log << c; return *this;}
  LogFile& operator<<(signed char c)                   {log << c; return *this;}
  LogFile& operator<<(const char *s)                   {log << s; return *this;}
  LogFile& operator<<(const unsigned char *s)          {log << s; return *this;}
  LogFile& operator<<(const signed char *s)            {log << s; return *this;}
  LogFile& operator<<(const void *p)                   {log << p; return *this;}
  LogFile& operator<<(int n)                           {log << n; return *this;}
  LogFile& operator<<(unsigned int n)                  {log << n; return *this;}
  // LogFile& operator<<(long n)                       {log << n; return *this;}
  // LogFile& operator<<(unsigned long n)              {log << n; return *this;}
  // LogFile& operator<<(short n)                      {log << n; return *this;}
  // LogFile& operator<<(unsigned short n)             {log << n; return *this;}
  LogFile& operator<<(double d)                        {log << d; return *this;}
  LogFile& operator<<(float d)                         {log << d; return *this;}
  LogFile& operator<<(bool b)                          {log << b; return *this;}
  LogFile& operator<<(std::string s)                   {log << s; return *this;}
  LogFile& operator<<(LogFile& (*f)(LogFile&))         {f(*this); return *this;}
  LogFile& operator>>(LogFile& (*f)(LogFile&))         {f(*this); return *this;}
  LogFile& flush()                                  {log.flush(); return *this;}
  void width(int n)                                              {log.width(n);}
  void fill(char n)                                               {log.fill(n);}
  void write(const Vector& v)                           {log << v << std::endl;}
  void write(const Matrix& m)                           {log << m << std::endl;}
};

// namespace Counters
// {
//   int c1;
//   int c2;
// };

void report(const double  d, const char* name="Noname",
            int total = 8, int decimal = 4);
void report(const Matrix& m, const char* name="Noname",
            int total = 8, int decimal = 4);
void report(const Vector& v,
            const char* name="Noname", int total = 8, int decimal = 4);

void add(Matrix& K, int row, int col, const Matrix& B1,
         const Matrix& C, const Matrix B2, double c1, double c0 = 0.);
void add(Vector& R, int row, const Matrix& BT,
         const Vector& V, double c1, double c0 = 0.);
void add2(Vector& R, int row, const Matrix& BT, const Vector& V,
          double c1, double c0 = 0.);

void add_BTCB(Matrix& K, int row, int col, const int* perm,
              const Matrix& B1, const Matrix& C, const Matrix B2,
              double c1, double c0 = 0.);
void add_BTv(Vector& R, int row, const int* perm,
             const Matrix& B, const Vector& v,
             double c1, double c0 = 0.);
void add_Bv(Vector& R, int row, const int* perm,
            const Matrix& B, const Vector& v,
            double c1, double c0 = 0.);

void spectralDecomposition(const Vector& s,
                           Vector& sP, Matrix& sV);
#endif  // SRC_MAIN_NEMESIS_DEBUG_H_
