#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "dixonTest.h"
/*
    Copyright (C) 2019 Thorsten Pohlert

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* define Fortran entry points, their names, nr of arguments*/
static const R_FortranMethodDef FortEntries[]  = {
  {"forpdixon", (DL_FUNC) &F77_SUB(forpdixon), 6},
  {"forqdixon", (DL_FUNC) &F77_SUB(forqdixon), 6},
  {"forddixon", (DL_FUNC) &F77_SUB(forddixon), 6},
  {NULL, NULL, 0}
};

// register
void attribute_visible R_init_dixonTest(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
