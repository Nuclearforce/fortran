C     PROGRAM SIGMA(INPUT,OUTPUT,TAPE1,TAPE2,TAPE3,TAPE6=OUTPUT)
      DIMENSION HLM(500),KVAR(500),IPMA(500)
      COMMON/MP/XP(500),ZP(500)
      COMMON/X/XJ(500,500),DEVD(500),DEVU(500),ELM(500),
     *SPIN(75),EN(75),IST(75),LEAD(2,500),VO(500),INP(75),
     *SJ2(80,5,5),SJ4(80,5,5),EM(75,75),VL(75),VP(75),T(4,75,75)
     *,IM(6),IP(6),CENV(4,14),SENV(14)
      COMMON/MVL/MV(14),ICO
      character*12 fil12,fil11,fil3
      mv(1)=2
      mv(2)=2
      mv(3)=2
      mv(4)=2
      mv(5)=4
      mv(6)=2
      mv(7)=4
      mv(8)=6
      mv(9)=2
      mv(10)=4
      mv(11)=4
      mv(12)=2
      mv(13)=4
      mv(14)=4
      DLTA=1.
      DO 33 I=1,14
      SENV(I)=0.
      DO 33 J=1,4
  33  CENV(J,I)=0.
      DO 100 I=1,500
      DO 100 J=1,500
  100 XJ(I,J)=0.
      DO 38 I=1,75
      IST(I)=0
      DO 38 J=1,75
 38   EM(I,J)=0.
      open(5,file='sigma.inp',status='old')
      READ(5,*)IPSW
      READ(5,*)NST
      read(5,811)fil12
      read(5,811)fil11
      read(5,811)fil3
  811 format(1a12)
      open(12,file=fil12,status='old')
      open(11,file=fil11,status='old')
      open(3,file=fil3,status='old')
      open(23,file='sigma.out',status='new')
      IF(NST.LE.0)GO TO 1
      NSTI=NST
      IF(NST.EQ.99)NSTI=75
      DO 2 I=1,NSTI
      J=I
      IF(NST.EQ.99)GO TO 2
      READ(5,*)J
  2   IST(J)=1
  1   READ(12,*)NMAX,MEMAX,INPO,INKO
      DO 3 I=1,NMAX
  3   READ(12,*)J,SPIN(J),EN(J)
      DO 4 I=1,MEMAX
  4   READ(12,*)J,LEAD(1,J),LEAD(2,J)
      DO 5 I=1,MEMAX
      READ(12,*)J,ELM(J)
  5   CONTINUE
  6   READ(12,*)J,I
      IF(J.EQ.0)GO TO 7
      READ(12,*)KH,IJ,D
      READ(12,*)(XJ(IPL,KH),IPL=1,MEMAX)
      READ(12,*)J,I
      READ(12,*)KH,IJ,DD
      DH=DD-D
      READ(12,*)(DEVD(IPL),IPL=1,MEMAX)
      XP(KH)=ELM(KH)-DD
      ZP(KH)=ELM(KH)-D
      DO 19 IPL=1,MEMAX
      XJ(IPL,KH)=(DEVD(IPL)-XJ(IPL,KH))/DH
  19  CONTINUE
      GO TO 6
  7   CONTINUE
      IF(NST.LT.0)GO TO 97
      READ(11,*)(DEVD(I),DEVU(I),I=1,MEMAX)
      DO 10 I=1,MEMAX
      IF(XJ(I,I).NE.0.)GO TO 17
       UUX=ABS(DEVD(I))
      UUY=DEVU(I)
      VV=AMAX1(UUX,UUY)
      DEVU(I)=VV
      DEVD(I)=-VV
      UUX=DEVD(I)
      UUY=DEVU(I)
      GO TO 18
  17  UUX=XP(I)*.6666667
      UUY=ZP(I)*.6666667
  18  CONTINUE
      VO(I)=-(DEVU(I)+DEVD(I))/2.
      XJ(I,I)=-2./UUX/UUY
      ZP(I)=1./XJ(I,I)
  10  CONTINUE
      DO 11 I=1,MEMAX
      DO 11 J=I,MEMAX
      IF(XJ(I,J).NE.0..AND.XJ(J,I).NE.0.)GO TO 14
      IF(XJ(I,J).EQ.0.)GO TO 15
      XJ(J,I)=XJ(I,J)
      GO TO 11
  15  XJ(I,J)=XJ(J,I)
      GO TO 11
  14  XJ(I,J)=.5*(XJ(I,J)+XJ(J,I))
      XJ(J,I)=XJ(I,J)
  11  CONTINUE
  97  CONTINUE
      DO 9 I=1,2
      DO 13 L=1,80
      DO 13 K=1,5
      IF(I.EQ.1)READ(3,*)(SJ2(L,K,J),J=1,5)
      IF(I.EQ.2)READ(3,*)(SJ4(L,K,J),J=1,5)
 13   CONTINUE
  9   CONTINUE
      DO 202 I=1,MEMAX
      HLM(I)=ELM(I)
 202  DEVD(I)=ELM(I)
      IF(NST.LT.0)GO TO 307
      IF(NST.NE.-1)CALL INVER(MEMAX,XJ)
      DO 303 I=1,MEMAX
      IF(XJ(I,I).GT.0.)GO TO 303
      XJ(I,I)=ZP(I)
  303 CONTINUE
  307  CONTINUE
      DO 203 IJA=1,NMAX
      DO 207 LL=1,MEMAX
  207 ELM(LL)=HLM(LL)
      ICO=0
      CALL INVR(IJA,NMAX,INPO,INKO,VAL)
      IF(NST.LT.0)GO TO 555
      DO 702 LL=1,MEMAX
 702  DEVD(LL)=HLM(LL)-VO(LL)
      IF(IPSW.NE.0)WRITE(23,888)
      DO 204 ICO=1,14
      IF(IST(IJA).NE.0)GO TO 8
      GO TO(8,8,8,8,8,204,204,204,8,8,8,204,204,204)ICO
  8   CONTINUE
      DO 703 LL=1,MEMAX
  703 ELM(LL)=DEVD(LL)
      CALL MINI(MEMAX,0,IJA,NMAX,INPO,INKO,DLTA,KVAR)
      DO 322 LL=1,MEMAX
  322 ELM(LL)=DEVD(LL)
      CALL MINI(MEMAX,1,IJA,NMAX,INPO,INKO,DLTA,KVAR)
      IF(IPSW.EQ.0)GO TO 204
      JLL=0
      DO 736 LL=INPO,INKO
      IF(KVAR(LL).EQ.0)GO TO 736
      JLL=JLL+1
      IPMA(JLL)=LL
 736  CONTINUE
      IF(ICO.EQ.1.OR.ICO.GE.9)IXCO=ICO
      IF(ICO.EQ.2)IXCO=8
      IF(ICO.GE.3.AND.ICO.LE.8)IXCO=ICO-1
      WRITE(23,735)IXCO
  735 FORMAT(1X,/5X,'CODE=',1I2)
      IF(JLL.NE.0)WRITE(23,*)(IPMA(LL),LL=1,JLL)
  204 CONTINUE
  555 CONTINUE
      WRITE(23,888)
      WRITE(23,700)
      DO 205 LLL=1,14
      EST=AMAX1(ABS(CENV(3,LLL)),ABS(CENV(4,LLL)))
      IF(CENV(3,LLL).GE.0.)CENV(3,LLL)=-EST
      IF(CENV(4,LLL).LE.0.)CENV(4,LLL)=EST
      IF(IST(IJA).NE.0)GO TO 205
      DO 217 LIL=1,3
      CENV(3,LIL+5)=0.
      CENV(4,LIL+5)=0.
      CENV(3,LIL+11)=0.
  217 CENV(4,LIL+11)=0.
 205  CONTINUE
 700  FORMAT(1X,125(1H*))
 888  FORMAT(1H1)
      WRITE(23,701)IJA,SPIN(IJA),EN(IJA)
      WRITE(23,700)
 701  FORMAT(5X,'INDEX=',1I2,5X,'SPIN=',1F4.1,5X,'ENERGY=',1F7.4,5X,/)
      WRITE(23,712)
      WRITE(23,718)CENV(1,1),CENV(3,1),CENV(4,1)
      DO 727 LL=1,3
      LI=2+LL
      LLL=2*LL-2
      WRITE(23,711)
      WRITE(23,713)LLL
      V=CENV(2,LI)
      IF(V.LT.0.)V=0.
      D=CENV(3,LI)
      U=CENV(4,LI)
      CALL SVAR(V,D,U)
 727  WRITE(23,719)CENV(1,LI),CENV(2,LI),CENV(3,LI),CENV(4,LI),V,D,U
     *,V/CENV(1,1),D/CENV(1,1),U/CENV(1,1)
      DO 728 LL=1,3
      LLL=2*LL-2
      LI=5+LL
      WRITE(23,711)
       WRITE(23,714)LLL
 728  WRITE(23,720)CENV(1,LI),CENV(2,LI),CENV(3,LI),CENV(4,LI)
      WRITE(23,711)
      WRITE(23,715)
      IF(CENV(1,3).LT.0.)CENV(1,3)=1.E-6
      WRITE(23,721)CENV(1,2),CENV(2,2),CENV(3,2),CENV(4,2),
     *FINT(1,CENV(1,1),CENV(1,3))
      DO 729 LL=1,3
      WRITE(23,711)
      LLL=2*LL-2
      LI=8+LL
      WRITE(23,716)LLL
      IF(CENV(1,2+LL).LT.0.)CENV(1,2+LL)=1.E-6
      IF(CENV(1,5+LL).LT.0.)CENV(1,5+LL)=1.E-6
      WRITE(23,721)CENV(1,LI),CENV(2,LI),CENV(3,LI),CENV(4,LI),
     *FINT(2,CENV(1,2+LL),CENV(1,5+LL))
 729  CONTINUE
      DO 730 LL=1,3
      LI=11+LL
      WRITE(23,711)
      V=CENV(2,LI)
      IF(V.LT.0.)V=0.
      D=CENV(3,LI)
      U=CENV(4,LI)
      CALL SVAR(V,D,U)
      WRITE(23,717)LL
 730  WRITE(23,719)CENV(1,LI),CENV(2,LI),CENV(3,LI),CENV(4,LI),V,D,U
  203 CONTINUE
      STOP
 718  FORMAT(5X,3(1F10.4,3X))
 719  FORMAT(5X,10(1F10.4,3X))
 720  FORMAT(5X,4(1F10.4,3X))
 721  FORMAT(5X,5(1F10.4,3X))
 711  FORMAT(1X,80(1H-))
 712  FORMAT(10X,'Q2',18X,'ERROR')
 713  FORMAT(10X,'Q4(',1I1,')',5X,'VARIANCE',15X,'ERROR',
     *14X,'SQRT(VAR)',11X,'ERROR',
     *11X,'SQRT(VAR)/Q2',11X,'ERROR')
 714  FORMAT(10X,'Q6(',1I1,')',5X,'SKEWNESS',15X,'ERROR')
 715  FORMAT(10X,'Q3COS(3D)',3X,'COS(3D)',14X,'ERROR',13X,'INT.Q3')
 716  FORMAT(10X,'Q5COS(3D)(',1I1,')',1X,'COS(3D)',13X,'ERROR',13X,
     *'INT.Q5')
 717  FORMAT(10X,'<COS2(3D)>(',1I1,')',1X,'VARIANCE',11X,'ERROR',12X,
     *'SQRT(VAR)',14X,'ERROR')
      END
      SUBROUTINE INVER(N,A)
      DIMENSION A(500,500),IK(500),JK(500)
      DO 100 K=1,N
      AM=0.
  21  CONTINUE
      DO 30 I=K,N
      DO 30 J=K,N
      IF(ABS(AM).GT.ABS(A(I,J)))GO TO 30
      AM=A(I,J)
      IK(K)=I
      JK(K)=J
  30  CONTINUE
      IF(AM.EQ.0.)GO TO 140
      I=IK(K)
      IF(I-K)21,51,43
  43  CONTINUE
      DO 50 J=1,N
      HOLD=A(K,J)
      A(K,J)=A(I,J)
  50  A(I,J)=-HOLD
  51  J=JK(K)
      IF(J-K)21,61,53
  53  CONTINUE
      DO 60 I=1,N
      HOLD=A(I,K)
      A(I,K)=A(I,J)
  60  A(I,J)=-HOLD
  61  CONTINUE
      DO 70 I=1,N
      IF(I.EQ.K)GO TO 70
      A(I,K)=-A(I,K)/AM
  70  CONTINUE
      DO 80 I=1,N
      DO 80 J=1,N
      IF(I.EQ.K)GO TO 80
      IF(J.EQ.K)GO TO 80
      A(I,J)=A(I,J)+A(I,K)*A(K,J)
  80  CONTINUE
      DO 90 J=1,N
      IF(J.EQ.K)GO TO 90
      A(K,J)=A(K,J)/AM
  90  CONTINUE
      A(K,K)=1./AM
 100  CONTINUE
      DO 130 L=1,N
      K=N-L+1
      J=IK(K)
      IF(K.GE.J)GO TO 111
      DO 110 I=1,N
      HOLD=A(I,K)
      A(I,K)=-A(I,J)
 110  A(I,J)=HOLD
 111  I=JK(K)
      IF(K.GE.I)GO TO 130
      DO 120 J=1,N
      HOLD=A(K,J)
      A(K,J)=-A(I,J)
 120  A(I,J)=HOLD
 130  CONTINUE
 140  RETURN
      END
      SUBROUTINE INVR(I,NMAX,INPO,INKO,VAL)
      DIMENSION IMA(14)
      COMMON/X/XJ(500,500),DEVD(500),DEVU(500),ELM(500),
     *SPIN(75),EN(75),IST(75),LEAD(2,500),VO(500),INP(75),
     *SJ2(80,5,5),SJ4(80,5,5),EM(75,75),VL(75),VP(75),T(4,75,75)
     *,IM(6),IP(6),CENV(4,14),SENV(14)
      COMMON/MVL/MV(14),ICO
      DO 18 J=INPO,INKO
      IJ=LEAD(1,J)
      IPL=LEAD(2,J)
      EM(IJ,IPL)=ELM(J)
  18  EM(IPL,IJ)=ELM(J)*(-1)**IFIX(SPIN(IJ)-SPIN(IPL))
      DO 17 J=1,NMAX
      VP(J)=EM(J,I)
  17  VL(J)=EM(I,J)
      NB=0
      DO 11 J=INPO,INKO
      IF(LEAD(1,J).NE.I.AND.LEAD(2,J).NE.I)GO TO 11
      NB=NB+1
      INP(NB)=J
  11  CONTINUE
      IF(ICO.NE.0)GO TO 6
      DO 1 L=1,14
      ICO=L
      IL=L
      DO 2 K=1,NMAX
      DO 2 M=1,NMAX
      DO 2 N=1,4
  2   T(N,K,M)=0.
      CALL INVAR(IL,I,NB,S,NMAX,INPO,INKO,SST)
       IF(L.EQ.12)S=S/SENV(6)
      IF(L.GE.13)S=6.*SST/SENV(8)-2.5*S/SENV(7)
      SENV(L)=S
  1   CENV(1,L)=S
      ICO=0
      GO TO 99
  6   CONTINUE
      IMA(1)=1
      GO TO(20,21,22,23,24,25,26,27,28,29,30,31,32,33)ICO
  20  ILE=1
      GO TO 88
  21  ILE=3
      IMA(2)=2
      IMA(3)=3
      GO TO 88
  22  ILE=2
      IMA(2)=3
      GO TO 88
  23  ILE=2
      IMA(2)=4
      GO TO 88
  24  ILE=2
      IMA(2)=5
      GO TO 88
  25  ILE=3
      IMA(2)=3
      IMA(3)=6
      GO TO 88
  26  ILE=3
      IMA(2)=4
      IMA(3)=7
      GO TO 88
  27  ILE=3
      IMA(2)=5
      IMA(3)=8
      GO TO 88
  28  ILE=4
      IMA(2)=3
      IMA(3)=6
      IMA(4)=9
      GO TO 88
  29  ILE=4
      IMA(2)=4
      IMA(3)=7
      IMA(4)=10
      GO TO 88
  30  ILE=4
      IMA(2)=5
      IMA(3)=8
      IMA(4)=11
      GO TO 88
  31  ILE=5
      IMA(2)=2
      IMA(3)=3
      IMA(4)=6
      IMA(5)=12
      GO TO 88
  32  ILE=6
      IMA(2)=2
      IMA(3)=3
      IMA(4)=7
      IMA(6)=8
      IMA(5)=13
      GO TO 88
  33  ILE=6
      IMA(2)=2
      IMA(3)=3
      IMA(4)=8
      IMA(6)=7
      IMA(5)=14
  88  CONTINUE
      DO 90 L=1,ILE
      IL=IMA(L)
      DO 89 K=1,NMAX
      DO 89 M=1,NMAX
      DO 89 N=1,4
  89  T(N,K,M)=0.
      ILI=IL
      CALL INVAR(ILI,I,NB,S,NMAX,INPO,INKO,SST)
      IF(IL.EQ.12)S=S/SENV(6)
      IF(IL.GE.13)S=6.*SST/SENV(8)-2.5*S/SENV(7)
  90  SENV(IL)=S
  99  CONTINUE
      IF(SENV(5).EQ.0.)SENV(5)=.00001
      IF(SENV(8).EQ.0.)SENV(8)=.00001
      IS=1
      IX=14
      IF(ICO.NE.0)IS=ICO
      IF(ICO.NE.0)IX=ICO
      DO 100 IGO=IS,IX
      VAL=SENV(1)
      IF(IGO.EQ.2)VAL=SENV(2)/FINT(1,SENV(1),SENV(3))
      IF(IGO.EQ.3)VAL=SENV(3)-SENV(1)*SENV(1)
      IF(IGO.EQ.4)VAL=SENV(4)-SENV(1)*SENV(1)
      IF(IGO.EQ.5)VAL=SENV(5)-SENV(1)*SENV(1)
      IF(IGO.EQ.6)VAL=SENV(6)-3.*SENV(3)*SENV(1)+2.*SENV(1)**3
      IF(IGO.EQ.7)VAL=SENV(7)-3.*SENV(4)*SENV(1)+2.*SENV(1)**3
      IF(IGO.EQ.8)VAL=SENV(8)-3.*SENV(5)*SENV(1)+2.*SENV(1)**3
      IF(IGO.EQ.9)VAL=SENV(9)/FINT(2,SENV(3),SENV(6))
      IF(IGO.EQ.10)VAL=SENV(10)/FINT(2,SENV(4),SENV(7))
      IF(IGO.EQ.11)VAL=SENV(11)/FINT(2,SENV(5),SENV(8))
      IF(IGO.EQ.12)VAL=SENV(12)-(SENV(2)/FINT(1,
     *SENV(1),SENV(3)))**2
      IF(IGO.EQ.13)VAL=SENV(13)-(SENV(2)/FINT(1,
     *SENV(1),SENV(3)))**2
      IF(IGO.EQ.14)VAL=SENV(14)-(SENV(2)/FINT(1,
     *SENV(1),SENV(3)))**2
      IF(ICO.EQ.0)CENV(2,IGO)=VAL
  100 CONTINUE
      RETURN
      END
      FUNCTION FINT(I,X,Y)
      IF(I.EQ.2)GO TO 2
      FINT=((SQRT(X)+SQRT(SQRT(Y)))/2.)**3
      RETURN
  2   Z=SQRT(SQRT(X))
      W=ALOG(Y)
      W=W/6.
      WW=EXP(W)
      Z=Z+WW
      Z=Z/2.
      Z=Z**5
      FINT=Z
      RETURN
      END
      SUBROUTINE INVAR(IN,I,NB,S,NMAX,INPO,INKO,SST)
      COMMON/X/XJ(500,500),DEVD(500),DEVU(500),ELM(500),
     *SPIN(75),EN(75),IST(75),LEAD(2,500),VO(500),INP(75),
     *SJ2(80,5,5),SJ4(80,5,5),EM(75,75),VL(75),VP(75),T(4,75,75)
     *,IM(6),IP(6),CENV(4,14),SENV(14)
      S=0.
      INTS=0
      IF(((-1)**IFIX(2.*SPIN(1)+.1)).NE.1)INTS=1
      Q=1./(2.*SPIN(I)+1.)
      GO TO(10,11,12,13,14,15,15,15,15,15,15,15,15,15)IN
 10   CONTINUE
      DO 1 J=1,NB
      L=INP(J)
  1   S=S+ELM(L)*ELM(L)
      S=S*Q
      RETURN
 11   CONTINUE
      DO 100 J=1,NB
      IL=INP(J)
      IL1=LEAD(1,IL)
      IF(ABS(SPIN(IL1)-SPIN(I)).GT.2.)GO TO 100
      IL=LEAD(2,IL)
      IF(ABS(SPIN(IL)-SPIN(I)).GT.2.)GO TO 100
      IV=IL1
      IF(IV.EQ.I)IV=IL
      DO 101 L=1,NMAX
      IF(ABS(SPIN(L)-SPIN(I)).GT.2.)GO TO 101
      IF((SPIN(L)+SPIN(I)).LT.2.)GO TO 101
      IF(EM(IV,L).EQ.0.)GO TO 101
      IT=IFIX(2.*SPIN(I)+.1)+1
      IU=IFIX(SPIN(L)-SPIN(I))+3
      ID=IFIX(SPIN(IV)-SPIN(I))+3
      T(2,L,IV)=EM(L,IV)*SJ2(IT,IU,ID)
      T(2,IV,L)=EM(IV,L)*SJ2(IT,IU,ID)
 101  CONTINUE
 100  CONTINUE
      IM(1)=2
      IP(1)=0
      CALL MULT(T,1,IM,VL,VP,S,IP,SPIN,I,NMAX)
      S=S*Q
      S=-S*4.183300133
      IF(INTS.EQ.1)S=-S
      RETURN
 12   CONTINUE
      IM(1)=1
      IM(2)=2
      IP(1)=1
      IP(2)=0
      CALL TJ(I,1,INPO,INKO)
      CALL MULT(T,2,IM,VL,VP,S,IP,SPIN,I,NMAX)
      S=5.*S*Q
      RETURN
  13  IM(1)=2
      IM(2)=2
      IP(1)=1
      IP(2)=0
      CALL TJ(I,2,INPO,INKO)
      CALL MULT(T,2,IM,VL,VP,S,IP,SPIN,I,NMAX)
      S=17.5*S*Q
      RETURN
 14   CONTINUE
      IM(1)=3
      IM(2)=1
      IP(1)=1
      IP(2)=0
      CALL TJ(I,3,INPO,INKO)
      CALL MULT(T,2,IM,VL,VP,S,IP,SPIN,I,NMAX)
      S=17.5*S*Q
      RETURN
  15  CONTINUE
      IN=IN-5
      DO 201 L=1,6
  201 IP(L)=0
      DO 16 J=1,NMAX
      DO 16 L=1,4
      DO 16 IV=1,NMAX
  16  T(L,J,IV)=0.
      GO TO(17,17,17,18,18,18,19,19,19)IN
  17  IM(1)=1
      IM(2)=4
      IM(3)=3
      J=(IN-1)*2
      CALL TEJ1(NMAX,I)
      CALL QJ(NMAX,I,J)
      CALL MULT(T,3,IM,VL,VP,S,IP,SPIN,I,NMAX)
      UU=17.5
      IF(J.EQ.0)UU=5.
      S=UU*Q*S
      RETURN
  18  J=(IN-4)*2
      DO 20 L=1,NMAX
  20  VP(L)=VL(L)
      CALL TEJ(NMAX,I,J)
      CALL QJ(NMAX,I,J)
      IM(1)=4
      IM(2)=3
      CALL MULT(T,2,IM,VL,VP,S,IP,SPIN,I,NMAX)
      UU=-73.20775232
      IF(J.EQ.0)UU=-20.91650066
      S=UU*S*Q
      IF(INTS.EQ.1)S=-S
      DO 21 L=1,NMAX
  21  VP(L)=EM(L,I)
      RETURN
  19  CONTINUE
      JO=-2
      NLX=1
      IF(IN.EQ.7)GO TO 55
      JO=0
      NLX=2
      IF(IN.EQ.9)GO TO 77
  55  CONTINUE
      DO 56 ILX=1,NLX
      J=JO+2*ILX
      IP(4)=2
      IP(2)=3
      IM(1)=2
      IM(2)=1
      IM(3)=3
      IM(4)=2
      CALL QJ(NMAX,I,2)
      IF(J.EQ.2)GO TO 31
      DO 44 L=1,NMAX
      DO 44 LL=1,NMAX
      T(3,L,LL)=0.
  44  T(1,L,LL)=0.
      CALL QJ(NMAX,I,J)
  31  CONTINUE
      CALL MULT(T,4,IM,VL,VP,AS,IP,SPIN,I,NMAX)
      AS=87.5*AS*Q
      IF(IN.EQ.7)S=AS
      IF(IN.EQ.8.AND.ILX.EQ.1)S=AS
      IF(IN.EQ.8.AND.ILX.EQ.2)SST=AS
  56  CONTINUE
      RETURN
  77   CONTINUE
      DO 57 ILX=1,NLX
      J=JO+2*ILX
      IP(1)=2
      IP(3)=2
      IM(1)=4
      IM(2)=3
      IM(3)=2
      CALL TEJ(NMAX,I,J)
      CALL QJ(NMAX,I,2)
      IF(J.EQ.2)GO TO 71
      DO 74 L=1,NMAX
      DO 74 LL=1,NMAX
      T(3,L,LL)=0.
  74  T(1,L,LL)=0.
      CALL QJ(NMAX,I,J)
  71  CONTINUE
      CALL MULT(T,3,IM,VL,VP,AS,IP,SPIN,I,NMAX)
      AS=AS*Q*87.5
      IF(ILX.EQ.1)S=AS
      IF(ILX.EQ.2)SST=AS
  57  CONTINUE
      RETURN
      END
      SUBROUTINE MULT(T,NM,IM,VL,VP,S,IP,SPIN,II,NMAX)
      DIMENSION T(4,75,75),IM(6),VL(75),VP(75),X(75),Y(75)
     *,IP(6),SPIN(75)
      COMMON/MVL/MV(14),ICO
      S=0.
      DO 1 I=1,NMAX
  1   Y(I)=VP(I)
      DO 2 I=1,NM
      DO 6 L=1,NMAX
  6   X(L)=0.
      L=IM(I)
      K=IP(I)
      DO 7 J=1,NMAX
      IF(IFIX(ABS(SPIN(J)-SPIN(II))).GT.MV(ICO))GO TO 7
      DO 3 IX=1,NMAX
  3   X(J)=T(L,J,IX)*PHASE(K,J,IX,SPIN,II)*Y(IX)+X(J)
  7   CONTINUE
      DO 4 J=1,NMAX
   4  Y(J)=X(J)
  2   CONTINUE
      DO 5 J=1,NMAX
  5   S=S+X(J)*VL(J)
      RETURN
      END
      FUNCTION PHASE(K,J,IX,SPIN,I)
      DIMENSION SPIN(75)
      PHASE=1.
      IF(K.EQ.0)RETURN
      GO TO (1,2,3)K
  1   PHASE=(-1)**IFIX(SPIN(I)+3.*SPIN(J)+.1)
      RETURN
  2   PHASE=(-1)**IFIX(SPIN(I)+SPIN(IX)+.1)
      RETURN
  3   PHASE=(-1)**IFIX(SPIN(J)+SPIN(IX)+.1)
      RETURN
      END
      SUBROUTINE TJ(I,K,INPO,INKO)
      COMMON/X/XJ(500,500),DEVD(500),DEVU(500),ELM(500),
     *SPIN(75),EN(75),IST(75),LEAD(2,500),VO(500),INP(75),
     *SJ2(80,5,5),SJ4(80,5,5),EM(75,75),VL(75),VP(75),T(4,75,75)
     *,IM(6),IP(6),CENV(4,14),SENV(14)
      S=SPIN(I)
      DO 2 L=INPO,INKO
      I1=LEAD(1,L)
      I2=LEAD(2,L)
      IF(EM(I1,I2).EQ.0.)GO TO 2
      IF(K.NE.1)GO TO 9
      DO 3 IV=1,2
      IF(IV.EQ.1)GO TO 4
      I1=LEAD(2,L)
      I2=LEAD(1,L)
  4   S1=SPIN(I1)
      S2=SPIN(I2)
      IF(S1.NE.S)GO TO 3
      T(1,I1,I2)=(-1)**IFIX(S+S2)*EM(I1,I2)
     */SQRT(5.*(2.*S+1.))
      T(2,I2,I1)=T(1,I1,I2)*(-1)**IFIX(S1-S2)
  3   CONTINUE
      GO TO 2
  9   CONTINUE
      IF(K.EQ.3)GO TO 19
      S1=SPIN(I1)
      S2=SPIN(I2)
      IF(ABS(S-S2).GT.2.)GO TO 2
      IF((S+S2).LT.2.)GO TO 2
      IN1=IFIX(2.*S2+1.1)
      IN2=IFIX(S-S2)+3
      IN3=IFIX(S1-S2)+3
      T(2,I1,I2)=EM(I1,I2)*SJ2(IN1,IN2,IN3)
      T(2,I2,I1)=EM(I2,I1)*SJ2(IN1,IN2,IN3)
 19   CONTINUE
      DO 20 IV=1,2
      IF(IV.EQ.1)GO TO 21
      I2=LEAD(1,L)
      I1=LEAD(2,L)
  21  CONTINUE
      S1=SPIN(I1)
      S2=SPIN(I2)
      IF(ABS(S-S1).GT.4.)GO TO 20
      IF(ABS(S-S2).GT.2.)GO TO 20
      IF((S+S1).LT.4.)GO TO 20
      IF((S+S2).LT.2.)GO TO 20
      IN1=IFIX(2.*S2+1.1)
      IN2=IFIX(S-S2)+3
      IN3=IFIX(S1-S2)+3
      T(3,I1,I2)=EM(I1,I2)*SJ4(IN1,IN2,IN3)
      T(1,I2,I1)=T(3,I1,I2)*(-1)**IFIX(S1-S2)
  20  CONTINUE
  2   CONTINUE
      RETURN
      END
      SUBROUTINE QJ(NMAX,I,J)
      COMMON/X/XJ(500,500),DEVD(500),DEVU(500),ELM(500),
     *SPIN(75),EN(75),IST(75),LEAD(2,500),VO(500),INP(75),
     *SJ2(80,5,5),SJ4(80,5,5),EM(75,75),VL(75),VP(75),T(4,75,75)
     *,IM(6),IP(6),CENV(4,14),SENV(14)
      S=SPIN(I)
      AJ=FLOAT(J)
      DO 1 K=1,NMAX
      S2=SPIN(K)
      IF(ABS(S-S2).GT.2.)GO TO 1
      IF((S+S2).LT.2.)GO TO 1
      IN1=IFIX(2.*S2+1.1)
      IN2=IFIX(S-S2)+3
      DO 5 L=1,NMAX
      IF(EM(K,L).EQ.0.)GO TO 5
      S1=SPIN(L)
      IF(ABS(S-S1).GT.AJ)GO TO 5
      IF((S+S1).LT.AJ)GO TO 5
      IN3=IFIX(S1-S2)+3
      IF(J.EQ.0)GO TO 2
      IF(J.EQ.2)GO TO 3
      T(3,K,L)=EM(K,L)*SJ4(IN1,IN2,IN3)
      T(1,L,K)=EM(L,K)*SJ4(IN1,IN2,IN3)
      GO TO 4
  2   T(3,K,L)=EM(K,L)*(-1)**IFIX(S+S2)
     */SQRT(5.*(2.*S+1.))
      T(1,L,K)=EM(L,K)*(-1)**IFIX(S+S2)
     */SQRT(5.*(2.*S+1.))
      GO TO 4
  3   T(3,K,L)=EM(K,L)*SJ2(IN1,IN2,IN3)
      T(1,L,K)=EM(L,K)*SJ2(IN1,IN2,IN3)
  4   CONTINUE
      IF(J.EQ.2)T(2,K,L)=T(3,K,L)
  5   CONTINUE
  1   CONTINUE
      RETURN
      END
      SUBROUTINE TEJ1(NMAX,I)
      COMMON/X/XJ(500,500),DEVD(500),DEVU(500),ELM(500),
     *SPIN(75),EN(75),IST(75),LEAD(2,500),VO(500),INP(75),
     *SJ2(80,5,5),SJ4(80,5,5),EM(75,75),VL(75),VP(75),T(4,75,75)
     *,IM(6),IP(6),CENV(4,14),SENV(14)
      DO 30 L=1,NMAX
      SJ=1./(2.*SPIN(L)+1.)
      DO 30 K=1,NMAX
      IF(SPIN(K).NE.SPIN(L))GO TO 30
      SS=0.
      DO 31 LU=1,NMAX
      IF(EM(LU,K).EQ.0..OR.EM(L,LU).EQ.0.)GO TO 31
      SS=SS+EM(L,LU)*EM(LU,K)*SJ*(-1)**IFIX(SPIN(LU)-SPIN(I))
  31  CONTINUE
      T(4,L,K)=SS
  30  CONTINUE
      RETURN
      END
      SUBROUTINE TEJ(NMAX,I,J)
      COMMON/X/XJ(500,500),DEVD(500),DEVU(500),ELM(500),
     *SPIN(75),EN(75),IST(75),LEAD(2,500),VO(500),INP(75),
     *SJ2(80,5,5),SJ4(80,5,5),EM(75,75),VL(75),VP(75),T(4,75,75)
     *,IM(6),IP(6),CENV(4,14),SENV(14)
      COMMON/MVL/MV(14),ICO
      S=SPIN(I)
      AJ=FLOAT(J)
      DO 5 L=1,NMAX
      S1=SPIN(L)
      IF(IFIX(ABS(S-S1)).GT.MV(ICO))GO TO 5
      IF(ABS(S-S1).GT.AJ)GO TO 5
      IF((S+S1).LT.AJ)GO TO 5
      DO 1 K=1,NMAX
      S2=SPIN(K)
      IF(ABS(S1-S2).GT.2.)GO TO 1
      IF((S1+S2).LT.2.)GO TO 1
      IF(ABS(S-S2).GT.2.)GO TO 1
      IF((S+S2).LT.2.)GO TO 1
      IN1=IFIX(2.*S2+1.1)
      IN3=IFIX(S1-S2)+3
      IN4=IFIX(S-S2)+3
      IF(J.EQ.0)GO TO 3
      IF(J.EQ.2)GO TO 4
      SJ=SJ4(IN1,IN3,IN4)
      GO TO 6
  4   SJ=SJ2(IN1,IN3,IN4)
      GO TO 6
  3   SJ=(-1)**IFIX(S+S2+.1)
     */SQRT(5.*(2.*S+1.))
  6   SS=0.
      DO 2 IV=1,NMAX
      IF(EM(L,IV).EQ.0..OR.EM(IV,K).EQ.0.)GO TO 2
      S3=SPIN(IV)
      IN2=IFIX(S3-S2)+3
      SSYM=SJ2(IN1,IN2,IN3)
      SS=SS+EM(L,IV)*EM(IV,K)*SSYM
  2   CONTINUE
      T(4,L,K)=SS*SJ
  1   CONTINUE
  5   CONTINUE
      RETURN
      END
      SUBROUTINE MINI(MEMAX,IDS,IJA,NMAX,INPO,INKO,DLTA,KVAR)
      DIMENSION KVAR(500),GRAD(500),GRAH(500),KFIX(500)
      COMMON/MP/XP(500),ZP(500)
      COMMON/X/XJ(500,500),DEVD(500),DEVU(500),ELM(500),
     *SPIN(75),EN(75),IST(75),LEAD(2,500),VO(500),INP(75),
     *SJ2(80,5,5),SJ4(80,5,5),EM(75,75),VL(75),VP(75),T(4,75,75)
     *,IM(6),IP(6),CENV(4,14),SENV(14)
      COMMON/MVL/MV(14),ICO
      NPTL=0
      IEF=0
      DO 888 I=1,MEMAX
      IF(KVAR(I).EQ.-2)KVAR(I)=1
      KFIX(I)=1
  888 CONTINUE
      IF(IDS.EQ.0)CHIL=9999.
      IF(IDS.EQ.1)CHIL=-9999.
  1   NPTL=NPTL+1
      DO 32 I=1,MEMAX
  32  DEVU(I)=ELM(I)
      CALL INVR(IJA,NMAX,INPO,INKO,CHISQ)
      IF(ABS(CHIL-CHISQ).LT..01)GO TO 999
       IF(NPTL.EQ.2)GO TO 22
      IF(IDS.EQ.0.AND.CHISQ.GE.CHIL)GO TO 77
      IF(IDS.EQ.1.AND.CHISQ.LE.CHIL)GO TO 77
  22  CHIL=CHISQ
      VALU=CHISQ-CENV(2,ICO)
      CENV(IDS+3,ICO)=VALU
      IF(IDS.NE.0)GO TO 2
      IF(NPTL.NE.1)GO TO 2
      DO 3 I=1,MEMAX
  3   KVAR(I)=1
  2   CONTINUE
      DO 4 I=INPO,INKO
      GRAD(I)=0.
      LMM=LEAD(1,I)
      LNN=LEAD(2,I)
      SSS=SPIN(LMM)
      TTT=SPIN(LNN)
      IF(IFIX(ABS(SPIN(IJA)-SSS)).GT.MV(ICO))KVAR(I)=0
      IF(IFIX(ABS(SPIN(IJA)-TTT)).GT.MV(ICO))KVAR(I)=0
      IF(KVAR(I).LE.0)GO TO 4
      XL=ELM(I)
      ELM(I)=ELM(I)+.01
      CALL INVR(IJA,NMAX,INPO,INKO,CHIS)
      GRAD(I)=100.*(CHIS-CHISQ)
      ELM(I)=XL
      IF(ABS(GRAD(I)).LT.1.E-6)KVAR(I)=0
      IF(NPTL.EQ.1)GRAH(I)=GRAD(I)
      IF(NPTL.EQ.1)GO TO 4
      IF((GRAD(I)*GRAH(I)).GT.0.)GO TO 4
      ELM(I)=(DEVD(I)*GRAD(I)-XL*GRAH(I))/(GRAD(I)-GRAH(I))
      GRAD(I)=0.
      KFIX(I)=0
      KVAR(I)=-2
  4   CONTINUE
  38  US=0.
      DO 10 I=1,MEMAX
      XP(I)=0.
      DO 11 II=1,MEMAX
      IF(IEF.EQ.1.AND.I.NE.II)GO TO 11
      XP(I)=XP(I)+XJ(I,II)*GRAD(II)
  11  CONTINUE
      US=US+XP(I)*GRAD(I)
  10  CONTINUE
      IF(US.GT.1.E-8)GO TO 37
      IF(IEF.EQ.1)GO TO 999
      IEF=1
      GO TO 38
  37  SS=1.
      IF(IDS.EQ.0)SS=-1.
      WSP=SS*SQRT(2.*DLTA/US)
      SS=0.
      DO 12 I=1,MEMAX
      IF(KFIX(I).EQ.0)GO TO 12
      US=WSP*XP(I)
      ELM(I)=DEVD(I)+US
  12  SS=SS+(ELM(I)-DEVU(I))**2
      SS=SQRT(SS)
      IF(SS.LT..005)GO TO 999
      GO TO 1
 999  CALL INVR(IJA,NMAX,INPO,INKO,CHISQ)
      VALU=CHISQ-CENV(2,ICO)
      CENV(IDS+3,ICO)=VALU
      RETURN
  77  CONTINUE
      IF(NPTL.LE.2.AND.IEF.EQ.0)GO TO 88
      RETURN
  88  IEF=1
      NPTL=0
      DO 44 I=1,MEMAX
  44  ELM(I)=DEVD(I)
      GO TO 1
      END
      SUBROUTINE SVAR(V,D,U)
      UG=V+U
      UD=V+D
      IF(UD.LT.0.)UD=1.E+20
      IF(UG.LT.0.)UG=1.E+20
      IF(V.LT.0.)V=1.E+10
      V=SQRT(V)
      D=SQRT(UD)-V
      U=SQRT(UG)-V
      RETURN
      END
