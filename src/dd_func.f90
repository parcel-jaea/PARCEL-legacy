module dd_func

contains
  ! c = a + b
  subroutine dd_add(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)
    implicit none
    real*8 a_hi,a_lo
    real*8 b_hi,b_lo
    real*8 c_hi,c_lo

    real*8 T1H,T1L
    real*8 X1H,X1L
    real*8 X2H,X2L
    real*8 p1,p2
    real*8 s1
    
    X1H  = a_hi
    X1L  = a_lo
    X2H  = b_hi
    X2L  = b_lo
    p1   = X1H + X2H
    p2   = p1 - X1H
    T1L  = (X1H - (p1 - p2)) + (X2H - p2)
    T1H  = p1
    T1L  = T1L + X1L + X2L
    s1   = T1H + T1L
    T1L  = T1L - (s1 - T1H)
    T1H  = s1
    c_hi = T1H
    c_lo = T1L

  end subroutine dd_add

  ! c = a - b
  subroutine dd_sub(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)
    implicit none
    real*8 a_hi,a_lo
    real*8 b_hi,b_lo
    real*8 c_hi,c_lo
    real*8 T1H,T1L
    real*8 X1H,X1L
    real*8 X2H,X2L
    real*8 p1,p2
    real*8 s1

    X1H  = a_hi
    X1L  = a_lo
    X2H  = b_hi
    X2L  = b_lo
    p1   = X1H - X2H
    p2   = p1 - X1H
    T1L  = (X1H - (p1 - p2)) - (X2H + p2)
    T1H  = p1
    T1L  = T1L + X1L - X2L
    s1   = T1H + T1L
    T1L  = T1L - (s1 - T1H)
    T1H  = s1
    c_hi = T1H
    c_lo = T1L

  end subroutine dd_sub

  ! c = a * b
  subroutine dd_mul(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)
    implicit none
    real*8 a_hi,a_lo
    real*8 b_hi,b_lo
    real*8 c_hi,c_lo

    real*8 T1H,T1L
    real*8 X1H,X1L
    real*8 X2H,X2L
    real*8 p1,p2
    real*8 r1,r2
    real*8 u1
    real*8 s1

    X1H  =  a_hi
    X1L  =  a_lo
    X2H  =  b_hi
    X2L  =  b_lo
    T1H  =  X1H * X2H
    u1   =  134217729.0d0*X1H
    p1   =  u1 - (u1-X1H)
    p2   =  X1H - p1
    u1   =  134217729.0d0*X2H
    r1   =  u1 - (u1-X2H)
    r2   =  X2H - r1
    T1L  =  ((p1*r1-T1H) + p1*r2 + p2*r1) + p2*r2
    T1L  =  T1L + (X1H * X2L + X1L * X2H)
    s1   =  T1H + T1L
    T1L  =  T1L - (s1 - T1H)
    T1H  =  s1
    c_hi = T1H
    c_lo = T1L

  end subroutine dd_mul

  ! c = a / b 
  subroutine dd_div(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo)
    implicit none
    real*8 a_hi,a_lo
    real*8 b_hi,b_lo
    real*8 c_hi,c_lo

    real*8 X1H,X1L
    real*8 X2H,X2L
    real*8 T1H,T1L
    real*8 p1,p2
    real*8 q1,q2
    real*8 r1,r2
    real*8 s1,s2
    real*8 u1
    
    X1H  = a_hi
    X1L  = a_lo
    X2H  = b_hi
    X2L  = b_lo
    q1   = X1H/X2H
    T1H  = X2H * q1
    u1   = 134217729.0d0*X2H
    p1   = u1 - (u1-X2H)
    p2   = X2H - p1
    u1   = 134217729.0d0*q1
    r1   = u1 - (u1-q1)
    r2   = q1 - r1
    T1L  = ((p1*r1-T1H) + p1*r2 + p2*r1) + p2*r2
    T1L  = T1L + (X2L * q1)
    s1   = T1H + T1L
    T1L  = T1L - (s1 - T1H)
    T1H  = s1
    s1   = X1H - T1H
    r1   = s1 - X1H
    s2   = (X1H - (s1 - r1)) - (T1H + r1)
    s2   = s2 - T1L
    s2   = s2 + X1L
    q2   = (s1 + s2)/X2H
    T1H  = q1 + q2
    T1L  = q2 - (T1H - q1)
    c_hi = T1H
    c_lo = T1L

  end subroutine dd_div

  ! b = sqrt(a)
  subroutine dd_sqrt(a_hi,a_lo,b_hi,b_lo)
    implicit none
    real*8 a_hi,a_lo
    real*8 b_hi,b_lo
    real*8 c_hi,c_lo

    real*8,parameter :: ONE=1.0d0
    real*8,parameter :: ZERO=0.0d0
    real*8,parameter :: HALF=0.5d0

    real*8 ABSXIH,SCALEH
    real*8 ABSXIL,SCALEL
    real*8 T1H,T2H
    real*8 T1L,T2L
    double precision s1,s2,q1,q2,r1,r2,p1,p2,t1,t2,u1,u2
    double precision x,ax,y1,y2,z1,z2

    if (a_hi.eq.0.0d0.and.a_lo.eq.0.0d0) then
       b_hi = ZERO
       b_lo = ZERO
    endif
    if (a_hi.lt.0.0d0.OR.(a_hi.eq.0.0d0.and.a_lo.lt.0.0d0)) then
       write(6,'(a)') 'SQRT CALLED WITH NEGATIVE ARGUMENT'
       b_hi = ZERO
       b_lo = ZERO
       return
    endif
    x = 1.0d0/sqrt(a_hi)
    ax = a_hi * x
    call dd_mul(x,ZERO,HALF,ZERO,y1,y2)
    ! y1 = x * HALF
    ! u1 = 134217729.0d0*x
    ! p1 = u1 - (u1-x)
    ! p2 = x - p1
    ! u1 = 134217729.0d0*HALF
    ! r1 = u1 - (u1-HALF)
    ! r2 = HALF - r1
    ! y2 = ((p1*r1-y1) + p1*r2 + p2*r1) + p2*r2
    ! s1 = y1 + y2
    ! y2 = y2 - (s1 - y1)
    ! y1 = s1


    call dd_mul(ax,ZERO,ax,ZERO,z1,z2)
    ! z1 = ax * ax
    ! u1 = 134217729.0d0*ax
    ! p1 = u1 - (u1-ax)
    ! p2 = ax - p1
    ! u1 = 134217729.0d0*ax
    ! r1 = u1 - (u1-ax)
    ! r2 = ax - r1
    ! z2 = ((p1*r1-z1) + p1*r2 + p2*r1) + p2*r2
    ! s1 = z1 + z2
    ! z2 = z2 - (s1 - z1)
    ! z1 = s1

    call dd_sub(a_hi,a_lo,z1,z2,T2H,T2L)
    ! p1 = a_hi - z1
    ! p2 = p1 - a_hi
    ! T2L = (a_hi - (p1 - p2)) - (z1 + p2)
    ! T2H = p1
    ! T2L = T2L + a_lo - z2
    ! s1 = T2H + T2L
    ! T2L = T2L - (s1 - T2H)
    ! T2H = s1

    call dd_mul(T2H,T2L,y1,y2,T1H,T1L)
    ! T1H = T2H * y1
    ! u1 = 134217729.0d0*T2H
    ! p1 = u1 - (u1-T2H)
    ! p2 = T2H - p1
    ! u1 = 134217729.0d0*y1
    ! r1 = u1 - (u1-y1)
    ! r2 = y1 - r1
    ! T1L = ((p1*r1-T1H) + p1*r2 + p2*r1) + p2*r2
    ! T1L = T1L + (T2H * y2 + T2L * y1)
    ! s1 = T1H + T1L
    ! T1L = T1L - (s1 - T1H)
    ! T1H = s1

    call dd_add(T1H,T1L,ax,ZERO, b_hi, b_lo)
    ! p1 = T1H + ax
    ! p2 = p1 - T1H
    ! q2 = (T1H - (p1 - p2)) + (ax - p2)
    ! q1 = p1
    ! q2 = q2 + T1L
    ! T1H = q1 + q2
    ! T1L = q2 - (T1H - q1)
    ! b_hi=T1H
    ! b_lo=T1L

  end subroutine dd_sqrt

end module dd_func

#if 0
! 4倍精度テスト
program main
  use dd_func
  implicit none
  integer i,ii
  integer,parameter ::n=32
  real*8 aa(n+2)
  real*8 aa_hi(n+2)
  real*8 aa_lo(n+2)
  real*8 aa0(n+2)
  real*8 aa1,aa2
  real*8 aa1_hi,aa1_lo
  real*8 aa2_hi,aa2_lo

  real*8 asqrt,bsqrt
  real*8 a_hi,a_lo
  real*8 b_hi,b_lo
  real*8 c_hi,c_lo

  real*8 a_true
  integer nn

  ! -----
  ! 漸化式 
  !   a(0)=1,a(1)=1/11,a(n+2)=34/11*a(n+1)-3/11*a(n)
  ! 一般解
  ! a(n)=1/11^n
  ! -----
  aa1=34.0d0/11.d0
  aa2=3.0d0/11.0d0
  call dd_div(34.0d0,0.0d0,11.d0,0.0d0,aa1_hi,aa1_lo)
  call dd_div(3.0d0 ,0.0d0,11.d0,0.0d0,aa2_hi,aa2_lo)

  aa(1)=1.0d0
  aa(2)=1.0d0/11.0d0

  aa_hi(1)=aa(1)
  aa_lo(1)=0.0d0
  call dd_div(1.0d0,0.0d0,11.d0,0.0d0,aa_hi(2),aa_lo(2))

  aa0(1)=1.0d0
  aa0(2)=1.0d0/11.0d0
  do i=3,n
     aa0(i)=aa0(i-1)/11.0d0
  enddo

  write(*,*)" True ,  double  ,  double - double"
  do i=3,n
     aa(i)=aa1*aa(i-1)-aa2*aa(i-2)
     call dd_mul(aa1_hi,aa1_lo,aa_hi(i-1),aa_lo(i-1),a_hi,a_lo)
     call dd_mul(aa2_hi,aa2_lo,aa_hi(i-2),aa_lo(i-2),b_hi,b_lo)
     call dd_sub(a_hi,a_lo,b_hi,b_lo,aa_hi(i),aa_lo(i))
     write(*,*)aa0(i),aa(i),aa_hi(i)
!     write(*,*)aa0(i),aa(i),aa_hi(i) + aa_lo(i)
  enddo

  ! -----
  ! 平方根
  ! -----
  nn=20
  a_true=1.1d0
  asqrt=a_true
  a_hi=asqrt
  a_lo=0.0d0
  do i=1,nn
     asqrt=sqrt(asqrt)
     call dd_sqrt(a_hi,a_lo,b_hi,b_lo) 
     a_hi=b_hi
     a_lo=b_lo
  enddo
  ! write(*,'(A20,E45.16)')       "double         ",asqrt
  ! write(*,'(A20,E45.16,E45.16)')"double - double",a_hi,a_lo

  ! -----
  ! ここで倍精度に直して4倍精度演算すると結果がずれるので4倍精度の平方根は正しい(下位ビットの情報が正しく計算できている)
  ! a_hi=asqrt
  ! a_lo=0.0d0
  do i=1,nn
     asqrt=asqrt*asqrt
     call dd_mul(a_hi,a_lo,a_hi,a_lo,b_hi,b_lo)
     a_hi=b_hi
     a_lo=b_lo
  enddo

  ! write(*,'(A20,E45.16)')       "True           ",a_true
  ! write(*,'(A20,E45.16)')       "double         ",asqrt
  ! write(*,'(A20,E45.16,E45.16)')"double - double",a_hi,a_lo

  write(*,*)"True           ",a_true
  write(*,*)"double         ",asqrt
  write(*,*)"double - double",a_hi
  write(*,*)" --------------------------- "

end program main
#endif
