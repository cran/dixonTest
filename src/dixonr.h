!  2019-05-14
!  T Pohlert
!  file named changed to dixonr.h
!  to comply with CRAN check policy

!  COMMON block for communication among Dixon-r routines.
!  Note that the very common variable names x, n, i, and j
!  are used in this COMMON block, for consistency with Dixon's notation;
!  therefore, subprograms that include this COMMON block must not
!  use those names for ordinary variables.

!  maximum numbers of quadrature points for r, x, and v integrations
      integer maxgl, maxfh, maxhh, maxnvec
      parameter (maxgl=16, maxfh=31, maxhh=17)
      parameter (maxnvec=maxfh*maxhh)

      common /dixonr/ xgl, wgl, t, u, x, w, z, c2, prefac, alpha_t,
     1     ngl, nfh, nhh, nvec, n, i, j

      double precision xgl(maxgl), wgl(maxgl), 
     1     t(maxnvec), u(maxnvec), x(maxnvec), w(maxnvec),
     2     z(maxnvec), c2(maxnvec), prefac, alpha_t

      integer ngl, nfh, nhh, nvec, n, i, j

