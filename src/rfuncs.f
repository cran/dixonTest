!  rfuncs.f
!  G. C. McBane (mcbaneg@gvsu.edu) January 2006
!  Journal of Statistical Software 16(3)
!  These programs evaluate the probability
!  density, cumulative distribution function, and critical values
!  for the r statistics of W. J. Dixon, The Annals of Mathematical 
!  Statistics, vol 22, pp. 68-78 (1951).


!  Function to return probability density function of one
!  of Dixon's r statistics.  The statistic is r_{j,i-1} for n
!  observations: r(i,j,n) = (x_n - x_{n-j})/(x_n - x_i),
!  where x_1 <= x_2 <= ... x_{n-1} <= x_n.
!  Values of j, i, and n must have been set in
!  a previous call to rinit or rreset.

      double precision function rdens(r)
      implicit none
      double precision r

      include 'dixonr.h'  !  common block containing i, j, etc.

      double precision fac1, f, v, g, c1, c3
      integer m

      double precision Phi  ! function for cumulative normal dist

!     Evaluate Eq. (11) of the text.  See subroutine
!     rinit for definitions of stored vectors.  m is a composite index
!     that runs over all quadrature points in t and u.

      fac1 = 1.0d0/sqrt(1+r*r)
      f = (1.0d0+r)*2*fac1*sqrt(1.0d0/3.0d0)
      rdens = 0.0d0
      do m = 1, nvec
         v = t(m)*fac1*sqrt(2.0d0)
         c1 = Phi(x(m) - v)
         c3 = Phi(x(m) - r*v)
!     The next line defining g has been replaced by the following four
!     lines to accomodate compilers that have difficulty with zero exponents
!        g = c1**(i-1)*(c3-c1)**(n-j-i-1)*(c2(m)-c3)**(j-1)
         g = 1.0d0
         if ((i-1).ne.0) g = g*c1**(i-1)
         if ((n-j-i-1).ne.0) g = g*(c3-c1)**(n-j-i-1)
         if ((j-1).ne.0) g = g*(c2(m)-c3)**(j-1)

         rdens = rdens + w(m)*v*exp(f*z(m))*g
      end do
      rdens = rdens*fac1*prefac*sqrt(4.0/3.0)
      return
      end

!     Return cumulative distribution function G(R) = prob(r<R) for
!     Dixon's r statistics.  Uses Gauss-Legendre quadrature.  Values of
!     i, j, n must have been set by previous call to rinit or rreset.

      double precision function rcdf(bigr)
      implicit none
      double precision bigr

      include 'dixonr.h'

      double precision rdens  !function
      double precision pr
      integer k
      
      rcdf = 0.0d0
      do k = 1, ngl
         pr = rdens(0.5*bigr*(xgl(k)+1))
         rcdf = rcdf+wgl(k)*pr
      end do
      rcdf = rcdf*bigr*0.5d0

      return
      end

!  Function to return critical values of R, that is, values of R such
!  that the cumulative distribution function G(R) = 1-alpha.
!  Uses root-finding routine zeroin() and auxilliary function
!  rcerr(R).  Values of i, j, n must have
!  been set by previous call to rinit or rreset.
      double precision function rcrit(alpha)
      implicit none
      double precision alpha

      double precision rcerr  !function to be passed to root-finder
      external rcerr

      double precision zeroin  ! root-finding function

      double precision errtol, a, b
      parameter(errtol = 1.0d-6, a = 0.0d0, b = 1.0d0)

      include 'dixonr.h'

      alpha_t = alpha

      rcrit = zeroin(a, b, rcerr, errtol)
      return
      end

!  Function called by zero-finding routine. Returns (1-alpha) - rcdf(R),
!  which will be zero when R is the appropriate critical value.

      double precision function rcerr(bigr)
      implicit none
      double precision bigr
      double precision rcdf !function

      include 'dixonr.h'

      rcerr = 1.0 - alpha_t - rcdf(bigr)

      return 
      end

!  Initializiation routine for Dixon r functions.  Must be
!  called once before any other functions in the package are used.
!  parameter prec can be 1 or 2 to select "low" or "high" precision;
!  it controls how many integration points are used in each
!  dimension.

      subroutine rinit(nt, it, jt, prec)
      implicit none
      integer nt, it, jt, prec

      include 'dixonr.h'

      double precision xfh(maxfh), wfh(maxfh), xhh(maxhh), whh(maxhh)
      integer m, l, k

      double precision Phi  ! function returning cumulative normal

      if (prec.eq.1) then
         ngl = maxgl-2
         nfh = maxfh-2
         nhh = maxhh-2
      else if (prec.eq.2) then
         ngl = maxgl
         nfh = maxfh
         nhh = maxhh
      end if

      nvec = nhh*nfh

      call hhquad(nhh,xhh,whh)  ! get half-range Hermite nodes and weights
      call fhquad(nfh, xfh, wfh) ! get full-range Hermite nodes and weights
      call glquad(ngl, xgl, wgl)  ! get Gauss-Legendre nodes and weights

!  Store vectors of r-independent quantities at the quadrature points in
!  t and u, for use by rdens().
      m = 0
      do l = 1, nhh             ! half-range index
         do k=1, nfh            ! full-range index
            m = m+1             ! composite index
            t(m) = xhh(l)       
            u(m) = xfh(k)
            x(m) = u(m)*sqrt(2.0d0/3.0d0)
            w(m) = whh(l)*wfh(k) ! combined weight
            z(m) = t(m)*u(m)     
            c2(m) = Phi(x(m))
         end do
      end do

      call rreset(nt,it,jt)

      return
      end

!  Store values of n, i, j, and 
!  compute and store normalization factor. 
      subroutine rreset(nt,it,jt)
      implicit none
      integer nt, it, jt

      double precision pi, pf
      integer k

      include 'dixonr.h'

!     Disabled 2019-05-14 by T Pohlert
!     to comply with CRAN policy and
!     R-ext manual       
!  Sanity checks

!      if ((nt.le.0).or.(it.le.0).or.(jt.le.0)) then
!         write(*,*) 'n, i, and j must all be greater than 0.'
!         write(*,*) 'n,i,j = ', nt, it, jt
!         stop
!      end if
!
!      if (nt.lt.(it+jt+1)) then
!         write(*,*) 'Inappropriate value of n.'
!         write(*,*) 'n,i,j = ', nt, it, jt
!         write(*,*) 'n should be larger than i+j.'
!         stop
!      end if

!  move requested values into common
      n = nt
      i = it
      j = jt

!  compute normalization factor (term in Dixon's eqn containing factorials,
!  plus three 1/sqrt(2*pi) terms from normal distributions)
      pi = acos(-1.0d0)
      prefac = sqrt(1.0d0/(2*pi)**3)
      do k = n, 1, -1
         pf = dble(k)
         if (k .le. i-1) pf = pf/k
         if (k.le. n-j-i-1) pf = pf/k
         if (k.le. j-1) pf = pf/k
         prefac = prefac*pf
      end do

      return
      end

