
c***********************************************************************
      SUBROUTINE LGNDR(PLM,MPLUS,LPLUS,THET)
c
c     Calculates Legendre polynomials Plm
c
c     mplus   number of m's >0
c     lplus   number of l's >0
c     thet    angle in degrees
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
 
      DIMENSION PLM(459)
c
c
c      write(6, *) "========== ",MPLUS, LPLUS, THET
c
      THETA=THET /57.295779
      Y=COS(THETA)
      Z=SIN(THETA)
      IX=0
      DO 100 M=1,MPLUS ! for MPLUS = 1, LPLUS = 16
      LX=M-1  ! LX = 0
      L2=0    ! L2 = 0
      P3=1.0  ! P3 = 1.0
      FL1=LX  ! FL1 = 0
      IF(LX.EQ.0) GO TO 41
      DO 40 LT=1,LX
      FL1=FL1+1.0
      P3=P3*FL1*Z/2.0
   40 CONTINUE
   41 P2=0.0  ! P2 = 0.0
      FL2=FL1+1.0  !FL2 = 1.0
      FL3=1.0    ! FL3 = 1.0
c================================= loop Lb
      DO 90 LT=1,LPLUS  ! loop Lb
      IX1=IX+LT 
      IF(L2-LX)50,70,60 ! if L2 < Lx -> 50; L2 == Lx -> 70; L2 > LX -> 60
   50 PLM(IX1)=0.0
      GO TO 75
   60 P3=(FL2*Y*P2-FL1*P1)/FL3
      FL1=FL1+1.0
      FL2=FL2+2.0
      FL3=FL3+1.0
   70 PLM(IX1)=P3

c      write(6, *) 'PLM, ',THETA*57.295779, IX1, PLM(IX1)

      P1=P2
      P2=P3
   75 L2=L2+1
   90 CONTINUE
c================================== end of Loop Lb
      IX=IX+LPLUS
  100 CONTINUE
      RETURN
      END
