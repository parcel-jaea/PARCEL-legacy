  ! -------------
  ! mul
  ! real*8 T1H,T1L
  ! real*8 X1H,X1L
  ! real*8 X2H,X2L
  ! real*8 p1,p2
  ! real*8 r1,r2
  ! real*8 u1
  ! real*8 s1

! #define ddmul_var@ T1H_ddmul_@,T1L_ddmul_@,X1H_ddmul_@,X1L_ddmul_@,X2H_ddmul_@,X2L_ddmul_@,p1_ddmul_@,p2_ddmul_@,r1_ddmul_@,r2_ddmul_@,u1_ddmul_@,s1_ddmul_@
#define ddmul_var1 T1H_ddmul_1,T1L_ddmul_1,X1H_ddmul_1,X1L_ddmul_1,X2H_ddmul_1,X2L_ddmul_1,p1_ddmul_1,p2_ddmul_1,r1_ddmul_1,r2_ddmul_1,u1_ddmul_1,s1_ddmul_1
#define ddmul_var2 T1H_ddmul_2,T1L_ddmul_2,X1H_ddmul_2,X1L_ddmul_2,X2H_ddmul_2,X2L_ddmul_2,p1_ddmul_2,p2_ddmul_2,r1_ddmul_2,r2_ddmul_2,u1_ddmul_2,s1_ddmul_2
#define ddmul_var3 T1H_ddmul_3,T1L_ddmul_3,X1H_ddmul_3,X1L_ddmul_3,X2H_ddmul_3,X2L_ddmul_3,p1_ddmul_3,p2_ddmul_3,r1_ddmul_3,r2_ddmul_3,u1_ddmul_3,s1_ddmul_3
#define ddmul_var4 T1H_ddmul_4,T1L_ddmul_4,X1H_ddmul_4,X1L_ddmul_4,X2H_ddmul_4,X2L_ddmul_4,p1_ddmul_4,p2_ddmul_4,r1_ddmul_4,r2_ddmul_4,u1_ddmul_4,s1_ddmul_4
#define ddmul_var5 T1H_ddmul_5,T1L_ddmul_5,X1H_ddmul_5,X1L_ddmul_5,X2H_ddmul_5,X2L_ddmul_5,p1_ddmul_5,p2_ddmul_5,r1_ddmul_5,r2_ddmul_5,u1_ddmul_5,s1_ddmul_5


  ! -------------
  ! X1H_ddmul_@  =  a_hi
  ! X1L_ddmul_@  =  a_lo
  ! X2H_ddmul_@  =  b_hi
  ! X2L_ddmul_@  =  b_lo
  ! T1H_ddmul_@  =  X1H_ddmul_@ * X2H_ddmul_@
  ! u1_ddmul_@   =  134217729.0d0*X1H_ddmul_@
  ! p1_ddmul_@   =  u1_ddmul_@ - (u1_ddmul_@-X1H_ddmul_@)
  ! p2_ddmul_@   =  X1H_ddmul_@ - p1_ddmul_@
  ! u1_ddmul_@   =  134217729.0d0*X2H_ddmul_@
  ! r1_ddmul_@   =  u1_ddmul_@ - (u1_ddmul_@-X2H_ddmul_@)
  ! r2_ddmul_@   =  X2H_ddmul_@ - r1_ddmul_@
  ! T1L_ddmul_@  =  ((p1_ddmul_@*r1_ddmul_@-T1H_ddmul_@) + p1_ddmul_@*r2_ddmul_@ + p2_ddmul_@*r1_ddmul_@) + p2_ddmul_@*r2_ddmul_@
  ! T1L_ddmul_@  =  T1L_ddmul_@ + (X1H_ddmul_@ * X2L_ddmul_@ + X1L_ddmul_@ * X2H_ddmul_@)
  ! s1_ddmul_@   =  T1H_ddmul_@ + T1L_ddmul_@
  ! T1L_ddmul_@  =  T1L_ddmul_@ - (s1_ddmul_@ - T1H_ddmul_@)
  ! T1H_ddmul_@  =  s1_ddmul_@
  ! c_hi = T1H_ddmul_@
  ! c_lo = T1L_ddmul_@
  ! #define ddmul_calc_@(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddmul_@=a_hi;X1L_ddmul_@=a_lo;X2H_ddmul_@=b_hi;X2L_ddmul_@=b_lo;T1H_ddmul_@=X1H_ddmul_@*X2H_ddmul_@;u1_ddmul_@=134217729.0d0*X1H_ddmul_@;p1_ddmul_@=u1_ddmul_@-(u1_ddmul_@-X1H_ddmul_@);p2_ddmul_@=X1H_ddmul_@-p1_ddmul_@;u1_ddmul_@=134217729.0d0*X2H_ddmul_@;r1_ddmul_@=u1_ddmul_@-(u1_ddmul_@-X2H_ddmul_@);r2_ddmul_@=X2H_ddmul_@-r1_ddmul_@;T1L_ddmul_@=((p1_ddmul_@*r1_ddmul_@-T1H_ddmul_@)+p1_ddmul_@*r2_ddmul_@+p2_ddmul_@*r1_ddmul_@)+p2_ddmul_@*r2_ddmul_@;T1L_ddmul_@=T1L_ddmul_@+(X1H_ddmul_@*X2L_ddmul_@+X1L_ddmul_@*X2H_ddmul_@);s1_ddmul_@=T1H_ddmul_@+T1L_ddmul_@;T1L_ddmul_@=T1L_ddmul_@-(s1_ddmul_@-T1H_ddmul_@);T1H_ddmul_@=s1_ddmul_@;c_hi=T1H_ddmul_@;c_lo=T1L_ddmul_@
  
#define ddmul_calc1(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddmul_1=a_hi;X1L_ddmul_1=a_lo;X2H_ddmul_1=b_hi;X2L_ddmul_1=b_lo;T1H_ddmul_1=X1H_ddmul_1*X2H_ddmul_1;u1_ddmul_1=134217729.0d0*X1H_ddmul_1;p1_ddmul_1=u1_ddmul_1-(u1_ddmul_1-X1H_ddmul_1);p2_ddmul_1=X1H_ddmul_1-p1_ddmul_1;u1_ddmul_1=134217729.0d0*X2H_ddmul_1;r1_ddmul_1=u1_ddmul_1-(u1_ddmul_1-X2H_ddmul_1);r2_ddmul_1=X2H_ddmul_1-r1_ddmul_1;T1L_ddmul_1=((p1_ddmul_1*r1_ddmul_1-T1H_ddmul_1)+p1_ddmul_1*r2_ddmul_1+p2_ddmul_1*r1_ddmul_1)+p2_ddmul_1*r2_ddmul_1;T1L_ddmul_1=T1L_ddmul_1+(X1H_ddmul_1*X2L_ddmul_1+X1L_ddmul_1*X2H_ddmul_1);s1_ddmul_1=T1H_ddmul_1+T1L_ddmul_1;T1L_ddmul_1=T1L_ddmul_1-(s1_ddmul_1-T1H_ddmul_1);T1H_ddmul_1=s1_ddmul_1;c_hi=T1H_ddmul_1;c_lo=T1L_ddmul_1

#define ddmul_calc2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddmul_2=a_hi;X1L_ddmul_2=a_lo;X2H_ddmul_2=b_hi;X2L_ddmul_2=b_lo;T1H_ddmul_2=X1H_ddmul_2*X2H_ddmul_2;u1_ddmul_2=134217729.0d0*X1H_ddmul_2;p1_ddmul_2=u1_ddmul_2-(u1_ddmul_2-X1H_ddmul_2);p2_ddmul_2=X1H_ddmul_2-p1_ddmul_2;u1_ddmul_2=134217729.0d0*X2H_ddmul_2;r1_ddmul_2=u1_ddmul_2-(u1_ddmul_2-X2H_ddmul_2);r2_ddmul_2=X2H_ddmul_2-r1_ddmul_2;T1L_ddmul_2=((p1_ddmul_2*r1_ddmul_2-T1H_ddmul_2)+p1_ddmul_2*r2_ddmul_2+p2_ddmul_2*r1_ddmul_2)+p2_ddmul_2*r2_ddmul_2;T1L_ddmul_2=T1L_ddmul_2+(X1H_ddmul_2*X2L_ddmul_2+X1L_ddmul_2*X2H_ddmul_2);s1_ddmul_2=T1H_ddmul_2+T1L_ddmul_2;T1L_ddmul_2=T1L_ddmul_2-(s1_ddmul_2-T1H_ddmul_2);T1H_ddmul_2=s1_ddmul_2;c_hi=T1H_ddmul_2;c_lo=T1L_ddmul_2

#define ddmul_calc3(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddmul_3=a_hi;X1L_ddmul_3=a_lo;X2H_ddmul_3=b_hi;X2L_ddmul_3=b_lo;T1H_ddmul_3=X1H_ddmul_3*X2H_ddmul_3;u1_ddmul_3=134217729.0d0*X1H_ddmul_3;p1_ddmul_3=u1_ddmul_3-(u1_ddmul_3-X1H_ddmul_3);p2_ddmul_3=X1H_ddmul_3-p1_ddmul_3;u1_ddmul_3=134217729.0d0*X2H_ddmul_3;r1_ddmul_3=u1_ddmul_3-(u1_ddmul_3-X2H_ddmul_3);r2_ddmul_3=X2H_ddmul_3-r1_ddmul_3;T1L_ddmul_3=((p1_ddmul_3*r1_ddmul_3-T1H_ddmul_3)+p1_ddmul_3*r2_ddmul_3+p2_ddmul_3*r1_ddmul_3)+p2_ddmul_3*r2_ddmul_3;T1L_ddmul_3=T1L_ddmul_3+(X1H_ddmul_3*X2L_ddmul_3+X1L_ddmul_3*X2H_ddmul_3);s1_ddmul_3=T1H_ddmul_3+T1L_ddmul_3;T1L_ddmul_3=T1L_ddmul_3-(s1_ddmul_3-T1H_ddmul_3);T1H_ddmul_3=s1_ddmul_3;c_hi=T1H_ddmul_3;c_lo=T1L_ddmul_3

#define ddmul_calc4(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddmul_4=a_hi;X1L_ddmul_4=a_lo;X2H_ddmul_4=b_hi;X2L_ddmul_4=b_lo;T1H_ddmul_4=X1H_ddmul_4*X2H_ddmul_4;u1_ddmul_4=134217729.0d0*X1H_ddmul_4;p1_ddmul_4=u1_ddmul_4-(u1_ddmul_4-X1H_ddmul_4);p2_ddmul_4=X1H_ddmul_4-p1_ddmul_4;u1_ddmul_4=134217729.0d0*X2H_ddmul_4;r1_ddmul_4=u1_ddmul_4-(u1_ddmul_4-X2H_ddmul_4);r2_ddmul_4=X2H_ddmul_4-r1_ddmul_4;T1L_ddmul_4=((p1_ddmul_4*r1_ddmul_4-T1H_ddmul_4)+p1_ddmul_4*r2_ddmul_4+p2_ddmul_4*r1_ddmul_4)+p2_ddmul_4*r2_ddmul_4;T1L_ddmul_4=T1L_ddmul_4+(X1H_ddmul_4*X2L_ddmul_4+X1L_ddmul_4*X2H_ddmul_4);s1_ddmul_4=T1H_ddmul_4+T1L_ddmul_4;T1L_ddmul_4=T1L_ddmul_4-(s1_ddmul_4-T1H_ddmul_4);T1H_ddmul_4=s1_ddmul_4;c_hi=T1H_ddmul_4;c_lo=T1L_ddmul_4

#define ddmul_calc5(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddmul_5=a_hi;X1L_ddmul_5=a_lo;X2H_ddmul_5=b_hi;X2L_ddmul_5=b_lo;T1H_ddmul_5=X1H_ddmul_5*X2H_ddmul_5;u1_ddmul_5=134217729.0d0*X1H_ddmul_5;p1_ddmul_5=u1_ddmul_5-(u1_ddmul_5-X1H_ddmul_5);p2_ddmul_5=X1H_ddmul_5-p1_ddmul_5;u1_ddmul_5=134217729.0d0*X2H_ddmul_5;r1_ddmul_5=u1_ddmul_5-(u1_ddmul_5-X2H_ddmul_5);r2_ddmul_5=X2H_ddmul_5-r1_ddmul_5;T1L_ddmul_5=((p1_ddmul_5*r1_ddmul_5-T1H_ddmul_5)+p1_ddmul_5*r2_ddmul_5+p2_ddmul_5*r1_ddmul_5)+p2_ddmul_5*r2_ddmul_5;T1L_ddmul_5=T1L_ddmul_5+(X1H_ddmul_5*X2L_ddmul_5+X1L_ddmul_5*X2H_ddmul_5);s1_ddmul_5=T1H_ddmul_5+T1L_ddmul_5;T1L_ddmul_5=T1L_ddmul_5-(s1_ddmul_5-T1H_ddmul_5);T1H_ddmul_5=s1_ddmul_5;c_hi=T1H_ddmul_5;c_lo=T1L_ddmul_5
