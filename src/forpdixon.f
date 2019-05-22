! Copyright (C) 2019 Thorsten Pohlert
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! A copy of the GNU General Public License is available at
! http://www.r-project.org/Licenses/

!     distribution function
      subroutine forpdixon(r, n, i, j, m, p)
      implicit none

      integer :: m
      double precision, dimension(m) :: p, r
      integer :: ii, i, j, n
      integer, parameter :: prec = 2
!     functions to be called
      double precision :: rcdf

!     initialise
      call rinit(n, i, j, prec)

!     calculate p
      do ii = 1, m
        p(ii) = rcdf(r(ii))
      end do

      end subroutine forpdixon

!    density functions
      subroutine forddixon(r, n, i, j, m, out)
      implicit none

      integer :: m
      double precision, dimension(m) :: out, r
      integer :: ii, i, j, n
      integer, parameter :: prec = 2
!     functions to be called
      double precision :: rdens

!     initialise
      call rinit(n, i, j, prec)

!     calculate p
      do ii = 1, m
        out(ii) = rdens(r(ii))
      end do

      end subroutine forddixon

!     quantile function
      subroutine forqdixon(p, n, i, j, m, r)
      implicit none

      integer :: m
      double precision, dimension(m) :: p, r
      integer :: ii, i, j, n
      integer, parameter :: prec = 2
!     functions to be called
      double precision :: rcrit

!     initialise
      call rinit(n, i, j, prec)

!     calculate p
      do ii = 1, m
        r(ii) = rcrit(p(ii))
      end do

      end subroutine forqdixon
