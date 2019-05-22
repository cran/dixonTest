#include <R.h>
#include <Rinternals.h>
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

void F77_NAME(forpdixon)(double *r, int *n, int *i, int *j, int *m, double *p);

void F77_NAME(forqdixon)(double *p, int *n, int *i, int *j, int *m, double *r);

void F77_NAME(forddixon)(double *r, int *n, int *i, int *j, int *m, double *out);
