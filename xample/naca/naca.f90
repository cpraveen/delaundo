!******************************************************************
!GENERATES POINTS ON THE BOUNDARY OF NACA 0012 
!******************************************************************
PROGRAM nacaprofile
   IMPLICIT NONE
   REAL,ALLOCATABLE,DIMENSION(:)::x,y,xo,yo
   REAL::theta,PI
   INTEGER::i,N,M,NoInnerBoundaries
   PRINT*,'Number of points on the profile'
   READ*,N
   PRINT*,'Number of points on the farfield boundary'
   READ*,M

   NoInnerBoundaries=1
   ALLOCATE(x(N),y(N),xo(M),yo(M))
   PI=4.0*ATAN(1.0)
   x(1)=1.00892
   y(1)=0.0
   DO i=1,N-1   !Finds points on the profile
      theta=2.0*PI*(i-1)/(N-2)
      x(i+1)=0.5044*(cos(theta)+1.0)
      y(i+1)=prof(x(i+1))
      IF(theta<PI)y(i+1)=-y(i+1)
   ENDDO

   DO i=1,N
      x(i)=x(i)/1.00892
      y(i)=y(i)/1.00892
   ENDDO


   DO i=1,M  !Calculates farfield points at a radius of 10 chord
      theta=2.0*PI*(i-1)/M   
      xo(i)=20.0*COS(theta)+0.5
      yo(i)=20.0*SIN(theta)
   ENDDO

   print*,'Writing airfoil data into naca.pts'
   OPEN(1,FILE='naca.pts')
   WRITE(1,*)'NEWBND'
   WRITE(1,*)'NAMEBN'
   WRITE(1,*) 1
   WRITE(1,*)'NRBNDE'
   WRITE(1,*) N-1
   WRITE(1,*)'NFRSBN'
   WRITE(1,*) 1
   WRITE(1,*)'NLSTBN'
   WRITE(1,*) 1
   WRITE(1,*)'ITYPBN'
   WRITE(1,*) 1
   WRITE(1,*)'BNDEXY'
   WRITE(1,5)x(1),y(1)
   WRITE(1,5)(x(i),y(i),i=3,N-1)
   WRITE(1,5)x(1),y(1)

   WRITE(1,*)'NEWBND'
   WRITE(1,*)'NAMEBN'
   WRITE(1,*) 2
   WRITE(1,*)'NRBNDE'
   WRITE(1,*) M+1
   WRITE(1,*)'NFRSBN'
   WRITE(1,*) 2
   WRITE(1,*)'NLSTBN'
   WRITE(1,*) 2
   WRITE(1,*)'ITYPBN'
   WRITE(1,*) 2
   WRITE(1,*)'BNDEXY'
   WRITE(1,5)(xo(i),yo(i),i=1,M)
   WRITE(1,5)xo(1),yo(1)
2  FORMAT(A10)
3  FORMAT(I5)
4  FORMAT(3I5)
6  FORMAT(2I5)
5  FORMAT(5X,2F15.6)
   CLOSE(1)

CONTAINS
   FUNCTION prof(x)
   IMPLICIT NONE
   REAL::x,C1,C2,C3,C4,C5
   REAL::prof

   C1=0.2969
   C2=-0.1260
   C3=-0.3516
   C4=0.2843
   C5=-0.1015
   prof=5.0*0.12*(C1*SQRT(x)+C2*x+C3*x**2+C4*x**3+C5*x**4)
   RETURN
   END FUNCTION prof

END PROGRAM nacaprofile

