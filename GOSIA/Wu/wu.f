      write(*,*)'A='
      read(*,*)a
      pi=2.*asin(1.)
   1  write(*,*)'Multipolarity='
      read(*,*)lam
      if(lam.eq.0)go to 99
      ie=1
      if(lam.le.6)go to 2
      ie=0
      lam=lam-6
   2  continue
      af=9./(lam+3)/(lam+3)
      aff=a**(2.*lam/3.)
      fff=.12**(2*lam)/pi
   3  write(*,*)'Me (>0 Wu to normal, <0 normal to Wu) ,spin (high)='
      read(*,*)xm,xs
      xs=2.*xs+1.
      if(xm.eq.0.)go to 1
      fac=fff*af*aff/4.
      if(ie.eq.0)fac=40.*fac/.12/.12/a**.66667
      write(*,*)'1Wu=',fac
      if(xm.lt.0.)go to 4
      s=abs(xm*xs*fac)
      s=sqrt(s)
      write(*,*)'Me=',s
      go to 3
  4   s=xm*xm/xs
      s=s/fac
      write(*,*)'Bl in Wu=',s
      go to 3
 99   continue
      stop
      end