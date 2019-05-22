!     2019-05-14
!     T Pohlert
!     write(,) and stop() expressions were disabled /
!     commented to comply with CRAN policy and
!     R-ext manual.
!
!  Fortran translation of Marsaglia's "little C function"
!  for the cumulative normal distribution
!  from Journal of Statistical Software 11(4)
!  Translation by G. C. McBane (mcbaneg@gvsu.edu) January 2006

      double precision function Phi(x)
      implicit none
      double precision x
      double precision s, t, b, q, i

      s = x
      t = 0.0d0
      b = x
      q = x*x
      i = 1.0d0

      do while (s .ne. t)
         t = s
         i = i+2.0d0
         b = b*q/i
         s = s+b
      end do
!  The obscure constant in the next line is ln(sqrt(2*pi)); it
!  is used instead of dividing the exponential by sqrt(2*pi)
      Phi = .5d0 + s*exp(-.5d0*q-.91893853320467274178d0)

      return
      end

!  Subroutine to return abscissas and weights for a Gauss-Hermite quadrature
!  on the half-infinite interval [0,infinity], for 29 or 31 quadrature points.
!  The tabular values  were generated
!  with the ORTHPOL package of W. Gautschi in double precision on an Intel
!  Pentium III.

      subroutine hhquad(n, x, w)
      implicit none
      integer n
      double precision x(n), w(n)
      integer i

      double precision xlow(15), wlow(15), xhigh(17), whigh(17)
      data xlow /
     *  0.2168694269885115D-01, 0.1126841969526419D+00,
     *  0.2704926203889978D+00, 0.4869022892342729D+00,
     *  0.7530435742305186D+00, 0.1060930872005834D+01,
     *  0.1404254809370576D+01, 0.1778646218520443D+01,
     *  0.2181707962845513D+01, 0.2613060672513553D+01,
     *  0.3074617939497498D+01, 0.3571407977467348D+01,
     *  0.4113735918455993D+01, 0.4723512894970656D+01,
     *  0.5460488773864854D+01 /
      data wlow /
     *  0.5544335429971267D-01, 0.1240277156956026D+00,
     *  0.1752909209537168D+00, 0.1914883325923960D+00,
     *  0.1634738094969453D+00, 0.1059376603780788D+00,
     *  0.5002704004358402D-01, 0.1644297800461658D-01,
     *  0.3573206793068966D-02, 0.4828969424775568D-03,
     *  0.3749090518569506D-04, 0.1493685973281760D-05,
     *  0.2552708581326519D-07, 0.1342178924114882D-09,
     *  0.9562291451849017D-13 /

      data xhigh /
     *  0.1809108332912565D-01, 0.9427561332967320D-01,
     *  0.2273676791408881D+00, 0.4116215432189120D+00,
     *  0.6404653042586674D+00, 0.9075572451603557D+00,
     *  0.1207442505307775D+01, 0.1535860903765161D+01,
     *  0.1889836898684654D+01, 0.2267684354439102D+01,
     *  0.2669032183793360D+01, 0.3094965841562238D+01,
     *  0.3548415607571656D+01, 0.4035068908990515D+01,
     *  0.4565576525100967D+01, 0.5161826321619212D+01,
     *  0.5882726076899502D+01 /
      data whigh /
     *  0.4629021825694398D-01, 0.1045285874343871D+00,
     *  0.1516928440386532D+00, 0.1752895232522240D+00,
     *  0.1651941135604728D+00, 0.1247589127315505D+00,
     *  0.7325004356830721D-01, 0.3228908736178809D-01,
     *  0.1029216269321119D-01, 0.2276385688421534D-02,
     *  0.3330379628636398D-03, 0.3037074862560045D-04,
     *  0.1594682391054344D-05, 0.4297674004833133D-07,
     *  0.4944542814585131D-09, 0.1723262832087000D-11,
     *  0.7782008254120001D-15 /

      if (n.eq.15) then
         do i = 1, n
            x(i) = xlow(i)
            w(i) = wlow(i)
         end do
      else if (n.eq.17) then
         do i = 1, n
            x(i) = xhigh(i)
            w(i) = whigh(i)
         end do
!      else
!         write(*,'(a,I3)') 'Illegal value of n in hhquad: n = ',! n
!         write(*,*) 'n must be 15 or 17.'
!         stop
      end if

      return
      end


!     Subroutine to return normal Gauss-Hermite weights and nodes for 29 or
!     31 quadrature points.
      subroutine fhquad(n,x,w)

      implicit none
      integer n
      double precision x(n), w(n)

      integer i
      double precision xlow(29), wlow(29), xhigh(31), whigh(31)
      data xlow /
     *  -.6728695198608847D+01, -.5998971289463819D+01,
     *  -.5389640521966748D+01, -.4841363651059163D+01,
     *  -.4331478293819151D+01, -.3848266792213620D+01,
     *  -.3384645141092214D+01, -.2935882504290126D+01,
     *  -.2498585691019404D+01, -.2070181076053428D+01,
     *  -.1648622913892315D+01, -.1232215755084753D+01,
     *  -.8194986812709121D+00, -.4091646363949289D+00,
     *  -.1105127916119648D-15, 0.4091646363949284D+00,
     *  0.8194986812709122D+00, 0.1232215755084752D+01,
     *  0.1648622913892316D+01, 0.2070181076053429D+01,
     *  0.2498585691019404D+01, 0.2935882504290126D+01,
     *  0.3384645141092216D+01, 0.3848266792213619D+01,
     *  0.4331478293819151D+01, 0.4841363651059162D+01,
     *  0.5389640521966747D+01, 0.5998971289463824D+01,
     *  0.6728695198608847D+01 /
      data wlow /
     *  0.1824460852767235D-19, 0.1534500444605336D-15,
     *  0.1390107271449587D-12, 0.3484130161308391D-10,
     *  0.3520312327600673D-08, 0.1749229129949956D-06,
     *  0.4823073497647784D-05, 0.7990920354521819D-04,
     *  0.8407925061402661D-03, 0.5845503545271516D-02,
     *  0.2763965559202364D-01, 0.9076884221557817D-01,
     *  0.2101426944492105D+00, 0.3464189390716706D+00,
     *  0.4089711746352303D+00, 0.3464189390716701D+00,
     *  0.2101426944492110D+00, 0.9076884221557793D-01,
     *  0.2763965559202356D-01, 0.5845503545271493D-02,
     *  0.8407925061402656D-03, 0.7990920354521827D-04,
     *  0.4823073497647759D-05, 0.1749229129949957D-06,
     *  0.3520312327600701D-08, 0.3484130161308462D-10,
     *  0.1390107271449590D-12, 0.1534500444605329D-15,
     *  0.1824460852767266D-19 /
      data xhigh /
     *  -.6995680123718540D+01, -.6275078704942860D+01,
     *  -.5673961444618588D+01, -.5133595577112374D+01,
     *  -.4631559506312859D+01, -.4156271755818146D+01,
     *  -.3700743403231469D+01, -.3260320732313541D+01,
     *  -.2831680453390201D+01, -.2412317705480420D+01,
     *  -.2000258548935638D+01, -.1593885860472139D+01,
     *  -.1191826998350046D+01, -.7928769769153088D+00,
     *  -.3959427364714228D+00, 0.4654056806286496D-15,
     *  0.3959427364714231D+00, 0.7928769769153078D+00,
     *  0.1191826998350046D+01, 0.1593885860472139D+01,
     *  0.2000258548935640D+01, 0.2412317705480417D+01,
     *  0.2831680453390206D+01, 0.3260320732313543D+01,
     *  0.3700743403231470D+01, 0.4156271755818143D+01,
     *  0.4631559506312856D+01, 0.5133595577112382D+01,
     *  0.5673961444618589D+01, 0.6275078704942855D+01,
     *  0.6995680123718541D+01 /
      data whigh /
     *  0.4618968394464040D-21, 0.5110609007927175D-17,
     *  0.5899556498753858D-14, 0.1860373521452180D-11,
     *  0.2352492003208647D-09, 0.1461198834491036D-07,
     *  0.5043712558939860D-06, 0.1049860275767554D-04,
     *  0.1395209039504714D-03, 0.1233683307306886D-02,
     *  0.7482799914035187D-02, 0.3184723073130041D-01,
     *  0.9671794816087061D-01, 0.2121327886687648D+00,
     *  0.3387726578941078D+00, 0.3957785560986103D+00,
     *  0.3387726578941073D+00, 0.2121327886687648D+00,
     *  0.9671794816087054D-01, 0.3184723073130041D-01,
     *  0.7482799914035165D-02, 0.1233683307306884D-02,
     *  0.1395209039504693D-03, 0.1049860275767565D-04,
     *  0.5043712558939813D-06, 0.1461198834491069D-07,
     *  0.2352492003208638D-09, 0.1860373521452120D-11,
     *  0.5899556498753877D-14, 0.5110609007927342D-17,
     *  0.4618968394464017D-21 /

      if (n.eq.29) then
         do i = 1, n
            x(i) = xlow(i)
            w(i) = wlow(i)
         end do
      else if (n.eq.31) then
         do i = 1, n
            x(i) = xhigh(i)
            w(i) = whigh(i)
         end do
!      else
!         write(*,'(a,I3)') 'Illegal value of n in fhquad: n = ', n
!         write(*,*) 'n must be 29 or 31.'
!         stop
      end if

      return
      end

!     Subroutine to return normal Gauss-Legendre weights and nodes on [-1,1] for
!     14 or 16 quadrature points.
      subroutine glquad(n,x,w)

      implicit none
      integer n
      double precision x(n), w(n)

      integer i
      double precision xlow(14), wlow(14), xhigh(16), whigh(16)
      data xlow /
     *  -.9862838086968126D+00, -.9284348836635737D+00,
     *  -.8272013150697650D+00, -.6872929048116854D+00,
     *  -.5152486363581541D+00, -.3191123689278898D+00,
     *  -.1080549487073434D+00, 0.1080549487073437D+00,
     *  0.3191123689278899D+00, 0.5152486363581542D+00,
     *  0.6872929048116855D+00, 0.8272013150697648D+00,
     *  0.9284348836635735D+00, 0.9862838086968123D+00/
      data wlow /
     *  0.3511946033175133D-01, 0.8015808715976042D-01,
     *  0.1215185706879034D+00, 0.1572031671581938D+00,
     *  0.1855383974779380D+00, 0.2051984637212957D+00,
     *  0.2152638534631579D+00, 0.2152638534631577D+00,
     *  0.2051984637212952D+00, 0.1855383974779380D+00,
     *  0.1572031671581937D+00, 0.1215185706879029D+00,
     *  0.8015808715976008D-01, 0.3511946033175180D-01/
      data xhigh /
     *  -.9894009349916500D+00, -.9445750230732329D+00,
     *  -.8656312023878321D+00, -.7554044083550031D+00,
     *  -.6178762444026439D+00, -.4580167776572275D+00,
     *  -.2816035507792586D+00, -.9501250983763712D-01,
     *  0.9501250983763745D-01, 0.2816035507792591D+00,
     *  0.4580167776572271D+00, 0.6178762444026440D+00,
     *  0.7554044083550032D+00, 0.8656312023878318D+00,
     *  0.9445750230732323D+00, 0.9894009349916499D+00 /
      data whigh /
     *  0.2715245941175416D-01, 0.6225352393864755D-01,
     *  0.9515851168249304D-01, 0.1246289712555336D+00,
     *  0.1495959888165769D+00, 0.1691565193950030D+00,
     *  0.1826034150449239D+00, 0.1894506104550681D+00,
     *  0.1894506104550689D+00, 0.1826034150449223D+00,
     *  0.1691565193950029D+00, 0.1495959888165762D+00,
     *  0.1246289712555355D+00, 0.9515851168249270D-01,
     *  0.6225352393864785D-01, 0.2715245941175414D-01 /

      if (n.eq.14) then
         do i = 1, n
            x(i) = xlow(i)
            w(i) = wlow(i)
         end do
      else if (n.eq.16) then
         do i = 1, n
            x(i) = xhigh(i)
            w(i) = whigh(i)
         end do
!      else
!         write(*,'(a,I3)') 'Illegal value of n in glquad: n = ', n
!         write(*,*) 'n must be 14 or 16.'
!         stop
      end if

      return
      end


      double precision function zeroin(ax,bx,f,tol)
      double precision ax,bx,f,tol
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result (.ge.0.)
c
c  output..
c
c  zeroin abscissa approximating a zero of  f  in the interval ax,bx
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  this is checked, and an error message is printed if this is not
c  satisfied.   zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
c  the  relative machine precision defined as the smallest representable
c  number such that  1.+macheps .gt. 1.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice-hall, inc. (1973).
c
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs, d1mach
      eps = d1mach(4)
c   10 eps = d1mach(4)
      tol1 = eps+1.0d0
c
      a=ax
      b=bx
      fa=f(a)
      fb=f(b)
c     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
c         write(6,2500)
c2500     format(1x,'f(ax) and f(bx) do not have different signs,',
c     1             ' zeroin is aborting')
c         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (dabs(fc).ge.dabs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol
      xm = 0.5d0*(c-b)
      if ((dabs(xm).le.tol1).or.(fb.eq.0.0d0)) go to 150
c
c see if a bisection is forced
c
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60
c
c linear interpolation
c
      p=2.0d0*xm*s
      q=1.0d0-s
      go to 70
c
c inverse quadratic interpolation
c
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
   70 if (p.le.0.0d0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.(p.ge.
     *dabs(0.5d0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (dabs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0d0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
  140 fb=f(b)
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
      go to 30
  150 zeroin=b
      return
      end

      DOUBLE PRECISION FUNCTION D1MACH(I)
      INTEGER I
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  D1MACH( 5) = LOG10(B)
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J
      COMMON /D9MACH/ CRAY1
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      DOUBLE PRECISION DMACH(5)
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
C  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
C  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
C  MANY MACHINES YET.
C  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
C  ON THE NEXT LINE
      DATA SC/0/
C  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
C  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
C          mail netlib@research.bell-labs.com
C          send old1mach from blas
C  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS.
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
C
C     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532
     *       .AND. SMALL(1) .EQ. -448790528) THEN
*           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935
     *       .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943
     *       .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074
     *       .AND. SMALL(2) .EQ. 58688) THEN
*           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 10               CONTINUE
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 20               CONTINUE
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
*                  *** CRAY ***
                  CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0)
                  SMALL(2) = 0
                  CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215)
                  CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214)
                  CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0)
                  RIGHT(2) = 0
                  CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0)
                  DIVER(2) = 0
                  CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215)
                  CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388)
               ELSE
c                  WRITE(*,9000)
c                  STOP 779
                  END IF
            ELSE
c               WRITE(*,9000)
c               STOP 779
               END IF
            END IF
         SC = 987
         END IF
*    SANITY CHECK
c      IF (DMACH(4) .GE. 1.0D0) STOP 778
c      IF (I .LT. 1 .OR. I .GT. 5) THEN
c         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
c         STOP
c         END IF
      D1MACH = DMACH(I)
      RETURN
c 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/
c     *' appropriate for your machine.')
* /* Standard C source for D1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*double d1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return DBL_MIN;
*	  case 2: return DBL_MAX;
*	  case 3: return DBL_EPSILON/FLT_RADIX;
*	  case 4: return DBL_EPSILON;
*	  case 5: return log10((double)FLT_RADIX);
*	  }
*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
*	exit(1); return 0; /* some compilers demand return values */
*}
      END
      SUBROUTINE I1MCRY(A, A1, B, C, D)
**** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END

