      program viscosity
c--------------------------------------------------------------------
c        Density functional theory for a case of 7B (HS) antibody model with
c          2 sticky sites of same stickiness - Wertheim association theory.       
c        PY approximation is taken to calculate the reference system.
c        ksum is the number of all sticky sites.
c        Authors: B. H. Lee & S. Brudar (April 2024)
c----------------------------------------------------->clean version

      implicit real*8(a-h,o-z)
      dimension x(4),hn(50),hhn(50),pn(50)
      character q2
      logical check
      common/hs/ sigm, omega
      common/n0/ pi
      common/deltaij/ delta_aa,delta_ab,delta_bb

      pi = 3.141592654d0
      cb = 1.38064852d-23
      h_p = 6.62607004d-34
      an = 6.02214179d23
      boltz_cal=3.297623483d-24
      ksum=2
      cm=0.0015d0
      delrho=1.d-6

      sigm=20.d0
      omega=1.8d0
      cc1=0.01205d0
      dd1=0.3762d0


      open(11,file='viscosity.in')
      open(17,file='viscosity.dat')

 200  format(2f10.5)     

      write(17,*)'MAB conc./(mg/mL)     viscosity/cP'



      read(11,*)q2,tt,q2,v175,q2,pm

      
      epsaa_cal=(dlog(v175)+4.902d0)/0.9588d0
      epsaa=epsaa_cal*1000.d0/(boltz_cal*an)
      epsbb=epsaa
      bepsaa = epsaa/tt
      bepsbb = epsbb/tt

c

      bepsab=dsqrt(bepsaa*bepsbb)
c
c***** Calculating deBroglie wavelength bl in A
      bl = h_p / dsqrt (2.d0 * pi * cb * tt * pm * 1d-3 / an)
      bl = bl * 1.d10

c
c###############################  rho loop ############################################      
c****** rho [/A**3] ************

      nmr=idint(cm/delrho)

      do ir=1,nmr

      cm1=1.d-10+dfloat(ir-1)*delrho

      rho = cm1 * an * 1.d-27

      rho_m=0.99 * rho
      rho_p=1.01 * rho

      etahs = 7.d0 * pi * rho * sigm**3 / 6.d0

      etahs_m = 7.d0 * pi * rho_m * sigm**3 / 6.d0
      etahs_p = 7.d0 * pi * rho_p * sigm**3 / 6.d0
c
c**** Calculating the ideal part of free energy beta*Aid/N - baid
      baid = 7.d0 * (dlog ( rho * bl**3) - 1.d0)

      baid_m = 7.d0 * (dlog ( rho_m * bl**3) - 1.d0)
      baid_p = 7.d0 * (dlog ( rho_p * bl**3) - 1.d0)

c**** Calculating the HS part of the free energy beta*Ahs/N -bahs
      bahs = 7.d0 * (etahs * (4.d0 - 3.d0 * etahs))/(1.d0 - etahs)**2

      bahs_m = 7.d0 * 
     *         (etahs_m * (4.d0 - 3.d0 * etahs_m))/(1.d0 - etahs_m)**2
      bahs_p = 7.d0 *
     *         (etahs_p * (4.d0 - 3.d0 * etahs_p))/(1.d0 - etahs_p)**2


c**** Wertheim
      call delta(rho,etahs,bepsaa,delta_ij)
      delta_aa=delta_ij
      call delta(rho,etahs,bepsab,delta_ij)
      delta_ab=delta_ij

c
c***** Calculating the fraction of sites ij that are not bonded - mass action law-


      x(1)=(-1.d0+dsqrt(1.d0+4.d0*dfloat(ksum)*rho*delta_aa))/
     /     (2.d0*dfloat(ksum)*rho*delta_aa)
      x(2)=x(1)


c*********************************************************************************
c**********************************************************************************


c
c***** Calculating the associative part of free energy beta*Aass/N - baass
      gsigma=(2.d0+etahs)/(2.d0*(1.d0-etahs)**2)
      baass=(dlog(x(1)*x(2)) - 
     *     0.5d0*(x(1)+x(2)+2.d0) + 4.d0*0.5d0) -
     *     6.d0*(dlog(rho*sigm**3*gsigma)-1.d0)


c
c**** 0.99 rho

      call delta(rho_m,etahs_m,bepsaa,delta_ij)
      delta_aa=delta_ij
      call delta(rho_m,etahs_m,bepsab,delta_ij)
      delta_ab=delta_ij

      x(1)=(-1.d0+dsqrt(1.d0+4.d0*dfloat(ksum)*rho_m*delta_aa))/
     /     (2.d0*dfloat(ksum)*rho_m*delta_aa)
      x(2)=x(1)


      gsigma_m=(2.d0+etahs_m)/(2.d0*(1.d0-etahs_m)**2)
      baass_m=(dlog(x(1)*x(2)) - 
     *     0.5d0*(x(1)+x(2)+2.d0) + 4.d0*0.5d0) -
     *     6.d0*(dlog(rho_m*sigm**3*gsigma_m)-1.d0)


c**** 1.01 rho

      call delta(rho_p,etahs_p,bepsaa,delta_ij)
      delta_aa=delta_ij
      call delta(rho_p,etahs_p,bepsab,delta_ij)
      delta_ab=delta_ij

      x(1)=(-1.d0+dsqrt(1.d0+4.d0*dfloat(ksum)*rho_p*delta_aa))/
     /     (2.d0*dfloat(ksum)*rho_p*delta_aa)
      x(2)=x(1)

      gsigma_p=(2.d0+etahs_p)/(2.d0*(1.d0-etahs_p)**2)
      baass_p=(dlog(x(1)*x(2)) - 
     *     0.5d0*(x(1)+x(2)+2.d0) + 4.d0*0.5d0) -
     *     6.d0*(dlog(rho_p*sigm**3*gsigma_p)-1.d0)



c
c***** Free energy of the system beta*A/N - bfree

       bfree = baid + bahs + baass

       bfree_m = baid_m + bahs_m + baass_m
       bfree_p = baid_p + bahs_p + baass_p

c***** Chemical potential (bchempot = beta*mu/N) and 
c***** osmotic pressure (osmp = beta*pressure /A**3)

       dif_m = (bfree*rho - bfree_m*rho_m)/(0.01 * rho)
       dif_p = (bfree_p*rho_p - bfree*rho)/(0.01 * rho)

       bchempot = (dif_m + dif_p)/2.d0
       bosmp = rho*bchempot - bfree*rho
        

c
c***** The probability distribution of n-size clusters hn(n)
       hn(1)=x(1)*x(2)
       hn(2)=(1.d0-x(1))*x(2)+
     +       (1.d0-x(2))*x(1) 
       hn(3)=(1.d0-x(1))*(1.d0-x(2))

c
c***** The viscosity calculation f(n) = cc1 * dfloat(n)**dd1
c***** The calculation is only valid for monoclonal antibodies : epsA=epsB=epsAB
c***** pn(j) is weight fraction distribution
c***** hhn(i) is n-distribution (probability that a cluster is a n-mere)

      avern=1.d0/x(1)
      relv=0.d0
      gammai=(rho/an/1.d-27)*pm
      do i=1,50
       hhn(i)=x(1)*(1.d0-x(1))**(dfloat(i)-1.d0)
       pn(i)=dfloat(i)*x(1)*hhn(i)
       relv=relv+(cc1*dfloat(i)**dd1)*pn(i)*gammai
      enddo


      do i=1,50
       hhn(i)=x(1)*(1.d0-x(1))**(dfloat(i)-1.d0)
       pn(i)=dfloat(i)*x(1)*hhn(i)
       relvn=(cc1*dfloat(i)**dd1)*pn(i)*gammai/relv
      enddo

       relvv=dexp(relv)
       write(17,200)gammai,relvv*1.0005d0

        

      enddo
c ############################ END rho loop #######################################
      close (12)

      end


c**********************************************************************************

      subroutine delta(rho,etahs,beps,delta_ij)
      implicit real*8(a-h,o-z)
      common/hs/ sigm, omega
      common/n0/ pi

      dd=sigm/2.d0

c***** Integral evaluation: fint = integral (f)       
       fint_an_p=4.d0*pi*(dexp(beps)-1.d0)/(6.d0*sigm**2)
       fint_int=(omega+sigm)**2 * (2.d0*omega-sigm) * 
     *             ((omega+sigm)**2/2.d0 - sigm**2/2.d0) +
     + ((omega+sigm)**2 - 2.d0*(omega+sigm)*(2.d0*omega-sigm))*
     *             ((omega+sigm)**3/3.d0 - sigm**3/3.d0) +
     + ((2.d0*omega-sigm)-2.d0*(omega+sigm))*
     *             ((omega+sigm)**4/4.d0 - sigm**4/4.d0) +
     *             ((omega+sigm)**5/5.d0 - sigm**5/5.d0) 
       fint=fint_int*fint_an_p
c
c***** Calculation of g(sigma) for HS in PY approximation
       gsigma=(2.d0+etahs)/(2.d0*(1.d0-etahs)**2)

c
c***** Calculation of Delta: Delta = g(sigma)* integral (f)
       delta_ij = gsigma * fint

       return
       end


c___________________________________________________________________________ 
      SUBROUTINE newt(x,n,rhot,check)
      implicit real*8(a-h,o-z)
      INTEGER n,nn,NP,MAXITS
      LOGICAL check
c      REAL x(n),fvec,TOLF,TOLMIN,TOLX,STPMX
c      REAL fvec,TOLF,TOLMIN,TOLX,STPMX
      PARAMETER (NP=40,MAXITS=200,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-7,
     *STPMX=100.)
      dimension x(n),fjac(NP,NP),g(NP),p(NP),xold(NP)
      COMMON /newtv/ fvec(NP),nn
c      common/deltaij/ delta_aa,delta_ab,delta_ac,delta_af,
c     *                delta_bb,delta_bc,delta_bf,
c     *                delta_cc,delta_cf,delta_ff
      common/deltaij/ delta_aa,delta_ab,delta_bb

      SAVE /newtv/
CU    USES fdjac,fmin,lnsrch,lubksb,ludcmp  
      INTEGER i,its,j,indx(NP)
c      REAL d,den,f,fold,stpmax,ssum,temp,test,fjac(NP,NP),g(NP),p(NP),
c     *xold(NP),fmin
      EXTERNAL fmin
      nn=n
c      write(*,*)'newt'
      f=fmin(x,rhot)
c      write(*,*)'f',f
      test=0.
      do 11 i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
11    continue
      if(test.lt..01*TOLF)then
        check=.false.
        return
      endif
      ssum=0.
      do 12 i=1,n
        ssum=ssum+x(i)**2
12    continue
      stpmax=STPMX*max(sqrt(ssum),float(n))
      do 21 its=1,MAXITS
        call fdjac(n,x,fvec,NP,rhot,fjac)
        do 14 i=1,n
          ssum=0.
          do 13 j=1,n
            ssum=ssum+fjac(j,i)*fvec(j)
13        continue
          g(i)=ssum
14      continue
        do 15 i=1,n
          xold(i)=x(i)
15      continue
        fold=f
        do 16 i=1,n
          p(i)=-fvec(i)
16      continue
c        write(*,*)'371'
        call ludcmp(fjac,n,NP,indx,d)
c        write(*,*)'373'
        call lubksb(fjac,n,NP,indx,p)
c        write(*,*)'375'
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,rhot,fmin)
c        write(*,*)'377'
        test=0.
        do 17 i=1,n
          if(abs(fvec(i)).gt.test)test=abs(fvec(i))
17      continue
        if(test.lt.TOLF)then
          check=.false.
          return
        endif
        if(check)then
          test=0.
          den=max(f,.5*n)
          do 18 i=1,n
            temp=abs(g(i))*max(abs(x(i)),1.)/den
            if(temp.gt.test)test=temp
18        continue
          if(test.lt.TOLMIN)then
            check=.true.
          else
            check=.false.
          endif
          return
        endif
        test=0.
        do 19 i=1,n
          temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
          if(temp.gt.test)test=temp
19      continue
        if(test.lt.TOLX)return
21    continue
c      pause 'MAXITS exceeded in newt'
      stop 'MAXITS exceeded in newt'
      END

      SUBROUTINE fdjac(n,x,fvec,np,rhot,df)
      implicit real*8(a-h,o-z)
      INTEGER n,np,NMAX
      PARAMETER (NMAX=4,EPS=1.e-4)
      dimension x(n),df(np,np),fvec(n),f(NMAX)
CU    USES funcv  
      INTEGER i,j
c      common/deltaij/ delta_aa,delta_ab,delta_ac,delta_af,
c     *                delta_bb,delta_bc,delta_bf,
c     *                delta_cc,delta_cf,delta_ff
      common/deltaij/ delta_aa,delta_ab,delta_bb

c      do 12 j=1,n  
c        temp=x(j)  
c        h=EPS*abs(temp)  
c        if(h.eq.0.)h=EPS  
c        x(j)=temp+h  
c        h=x(j)-temp  
c        call funcv(n,x,f)  
c        x(j)=temp  
c        do 11 i=1,n  
c          df(i,j)=(f(i)-fvec(i))/h  
c11      continue  
c12    continue  

      df(1,1)=1.d0+2.d0*rhot*delta_aa*x(1)+rhot*delta_ab*x(2)+
     +                  rhot*delta_ac*x(3)+rhot*delta_af*x(4)
      df(1,2)=rhot*delta_ab*x(1)
      df(1,3)=rhot*delta_ac*x(1)
      df(1,3)=rhot*delta_af*x(1)

      df(2,1)=rhot*delta_ab*x(2)
      df(2,2)=1.d0+rhot*delta_ab*x(1)+2.d0*rhot*delta_bb*x(2)+
     +                  rhot*delta_bc*x(3)+rhot*delta_bf*x(4)
      df(3,2)=rhot*delta_bc*x(2)
      df(4,2)=rhot*delta_bf*x(2)

      df(3,1)=rhot*delta_ac*x(3)
      df(3,2)=rhot*delta_bc*x(3)
      df(3,3)=1.d0+rhot*delta_ac*x(1)+rhot*delta_bc*x(2)+
     +                  2.d0*rhot*delta_cc*x(3)+rhot*delta_cf*x(4)
      df(3,4)=rhot*delta_cf*x(3)

      df(4,1)=rhot*delta_af*x(4)
      df(4,2)=rhot*delta_bf*x(4)
      df(4,3)=rhot*delta_cf*x(4)
      df(4,4)=1.d0+rhot*delta_af*x(1)+rhot*delta_bf*x(2)+
     +                  rhot*delta_cf*x(3)+2.d0*rhot*delta_ff*x(4)


      return
      END


      subroutine funcv(n,x,rhot,f)
      implicit real*8(a-h,o-z)
      dimension x(n),f(n)
c      common/deltaij/ delta_aa,delta_ab,delta_ac,delta_af,
c     *                delta_bb,delta_bc,delta_bf,
c     *                delta_cc,delta_cf,delta_ff
      common/deltaij/ delta_aa,delta_ab,delta_bb

      f(1)=x(1)+rhot*delta_aa*x(1)*x(1)+rhot*delta_ab*x(1)*x(2)+
     +          rhot*delta_ac*x(1)*x(3)+rhot*delta_af*x(1)*x(4)-1.d0

      f(2)=x(2)+rhot*delta_ab*x(1)*x(2)+rhot*delta_bb*x(2)*x(2)+
     +          rhot*delta_bc*x(2)*x(3)+rhot*delta_bf*x(2)*x(4)-1.d0

      f(3)=x(3)+rhot*delta_ac*x(1)*x(3)+rhot*delta_bc*x(2)*x(3)+
     +          rhot*delta_cc*x(3)*x(3)+rhot*delta_cf*x(3)*x(4)-1.d0

      f(4)=x(4)+rhot*delta_af*x(1)*x(4)+rhot*delta_bf*x(2)*x(4)+
     +          rhot*delta_cf*x(3)*x(4)+rhot*delta_ff*x(4)*x(4)-1.d0
c      write(*,*)'funcv',f(1),f(2),f(3),f(4)

      return
      end


      FUNCTION fmin(x,rhot)
      implicit real*8(a-h,o-z)
      INTEGER n,NP
      PARAMETER (NP=40)
      COMMON /newtv/ fvec(NP),n
      dimension x(n)
      SAVE /newtv/
CU    USES funcv  
      INTEGER i
c      write(*,*)'fmin'
      call funcv(n,x,rhot,fvec)
      ssum=0.
      do 11 i=1,n
        ssum=ssum+fvec(i)**2
11    continue
      fmin=0.5*ssum
      return
      END


      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,rhot,func)
      implicit real*8(a-h,o-z)
      INTEGER n
      LOGICAL check
      dimension x(n),g(n),p(n),xold(n)
      PARAMETER (ALF=1.e-4,TOLX=1.e-7)
      EXTERNAL func
CU    USES func  
      INTEGER i
c      write(*,*)'lnsrch'
      check=.false.
      ssum=0.
      do 11 i=1,n
        ssum=ssum+p(i)*p(i)
11    continue
      ssum=sqrt(ssum)
      if(ssum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/ssum
12      continue
      endif
      slope=0.
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.
c      write(*,*)'540'
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
c        write(*,*)'545',x(1),x(2),x(3),x(4),f
        f=func(x,rhot)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
c          write(*,*)'550'
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.)then
            tmplam=-slope/(2.*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.)then
              tmplam=-slope/(2.*b)
            else
              disc=b*b-3.*a*slope
c              if(disc.lt.0.) pause 'roundoff problem in lnsrch'
              if(disc.lt.0.) stop 'roundoff problem in lnsrch'
              tmplam=(-b+sqrt(disc))/(3.*a)
            endif
            if(tmplam.gt..5*alam)tmplam=.5*alam
          endif

c          write(*,*)'574'
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,.1*alam)
      goto 1
      END



      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      implicit real*8(a-h,o-z)
      DIMENSION A(NP,NP),INDX(N),B(N)
c      write(*,*)'lubksb'
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END



      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      implicit real*8(a-h,o-z)
      PARAMETER (NMAX=100,TINY=1.0E-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
c        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        IF (AAMAX.EQ.0.) stop 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
