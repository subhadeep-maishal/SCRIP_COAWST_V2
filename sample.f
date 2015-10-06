      Program check 
      real :: xx(4), px
      real :: yy(4), py 
      integer :: N

      xx(1) = 0.0 
      xx(2) = 5.0 
      xx(3) = 5.0
      xx(4) = 0.0 
      
      yy(1) = 0.0 
      yy(2) = 0.0 
      yy(3) = 5.0
      yy(4) = 5.0

      do i = 1,4
      write(20,*) xx(i),yy(i)
      end do 

      px=0.0
      py=5.0 
      N =4
      
      call PNPOLY(PX,PY,XX,YY,N,INOUT) 
      print*, INOUT                     
      End program check
        
      SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT)                            
      REAL X(200),Y(200),XX(N),YY(N)                                    
      LOGICAL MX,MY,NX,NY                                               
      INTEGER O                                                         
C      OUTPUT UNIT FOR PRINTED MESSAGES                                 
      DATA O/6/                                                         
      MAXDIM=200                                                        
      IF(N.LE.MAXDIM)GO TO 6                                            
      WRITE(O,7)                                                        
7     FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY.     
     1RESULTS INVALID')                                                 
      RETURN                                                            
6     DO 1 I=1,N                                                        
      X(I)=XX(I)-PX                                                     
1     Y(I)=YY(I)-PY                                                     
      INOUT=-1                                                          
      DO 2 I=1,N                                                        
      J=1+MOD(I,N)                                                      
      MX=X(I).GE.0.0                                                    
      NX=X(J).GE.0.0                                                    
      MY=Y(I).GE.0.0                                                    
      NY=Y(J).GE.0.0                                                    
      IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2       
      IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
      INOUT=-INOUT                                                      
      GO TO 2                                                           
3     IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5                       
4     INOUT=0                                                           
      RETURN                                                            
5     INOUT=-INOUT                                                      
2     CONTINUE                                                          
      RETURN                                                            
      END                             
