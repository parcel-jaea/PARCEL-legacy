! -------------
! div
!    real*8 X1H,X1L
!    real*8 X2H,X2L
!    real*8 T1H,T1L
!    real*8 p1,p2
!    real*8 q1,q2
!    real*8 r1,r2
!    real*8 s1,s2
!    real*8 u1

! #define dddiv_var@ X1H_dddiv_@,X1L_dddiv_@,X2H_dddiv_@,X2L_dddiv_@,T1H_dddiv_@,T1L_dddiv_@,p1_dddiv_@,p2_dddiv_@,q1_dddiv_@,q2_dddiv_@,r1_dddiv_@,r2_dddiv_@,s1_dddiv_@,s2_dddiv_@,u1_dddiv_@

#define dddiv_var1 X1H_dddiv_1,X1L_dddiv_1,X2H_dddiv_1,X2L_dddiv_1,T1H_dddiv_1,T1L_dddiv_1,p1_dddiv_1,p2_dddiv_1,q1_dddiv_1,q2_dddiv_1,r1_dddiv_1,r2_dddiv_1,s1_dddiv_1,s2_dddiv_1,u1_dddiv_1

#define dddiv_var2 X1H_dddiv_2,X1L_dddiv_2,X2H_dddiv_2,X2L_dddiv_2,T1H_dddiv_2,T1L_dddiv_2,p1_dddiv_2,p2_dddiv_2,q1_dddiv_2,q2_dddiv_2,r1_dddiv_2,r2_dddiv_2,s1_dddiv_2,s2_dddiv_2,u1_dddiv_2

#define dddiv_var3 X1H_dddiv_3,X1L_dddiv_3,X2H_dddiv_3,X2L_dddiv_3,T1H_dddiv_3,T1L_dddiv_3,p1_dddiv_3,p2_dddiv_3,q1_dddiv_3,q2_dddiv_3,r1_dddiv_3,r2_dddiv_3,s1_dddiv_3,s2_dddiv_3,u1_dddiv_3

#define dddiv_var4 X1H_dddiv_4,X1L_dddiv_4,X2H_dddiv_4,X2L_dddiv_4,T1H_dddiv_4,T1L_dddiv_4,p1_dddiv_4,p2_dddiv_4,q1_dddiv_4,q2_dddiv_4,r1_dddiv_4,r2_dddiv_4,s1_dddiv_4,s2_dddiv_4,u1_dddiv_4

#define dddiv_var5 X1H_dddiv_5,X1L_dddiv_5,X2H_dddiv_5,X2L_dddiv_5,T1H_dddiv_5,T1L_dddiv_5,p1_dddiv_5,p2_dddiv_5,q1_dddiv_5,q2_dddiv_5,r1_dddiv_5,r2_dddiv_5,s1_dddiv_5,s2_dddiv_5,u1_dddiv_5

! -------------
!     X1H  = a_hi
!     X1L  = a_lo
!     X2H  = b_hi
!     X2L  = b_lo
!     q1   = X1H/X2H
!     T1H  = X2H * q1
!     u1   = 134217729.0d0*X2H
!     p1   = u1 - (u1-X2H)
!     p2   = X2H - p1
!     u1   = 134217729.0d0*q1
!     r1   = u1 - (u1-q1)
!     r2   = q1 - r1
!     T1L  = ((p1*r1-T1H) + p1*r2 + p2*r1) + p2*r2
!     T1L  = T1L + (X2L * q1)
!     s1   = T1H + T1L
!     T1L  = T1L - (s1 - T1H)
!     T1H  = s1
!     s1   = X1H - T1H
!     r1   = s1 - X1H
!     s2   = (X1H - (s1 - r1)) - (T1H + r1)
!     s2   = s2 - T1L
!     s2   = s2 + X1L
!     q2   = (s1 + s2)/X2H
!     T1H  = q1 + q2
!     T1L  = q2 - (T1H - q1)
!     c_hi = T1H
!     c_lo = T1L

! #define dddiv_calc@(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_dddiv_@=a_hi;X1L_dddiv_@=a_lo;X2H_dddiv_@=b_hi;X2L_dddiv_@=b_lo;q1_dddiv_@=X1H_dddiv_@/X2H_dddiv_@;T1H_dddiv_@=X2H_dddiv_@*q1_dddiv_@;u1_dddiv_@=134217729.0d0*X2H_dddiv_@;p1_dddiv_@=u1_dddiv_@-(u1_dddiv_@-X2H_dddiv_@);p2_dddiv_@=X2H_dddiv_@-p1_dddiv_@;u1_dddiv_@=134217729.0d0*q1_dddiv_@;r1_dddiv_@=u1_dddiv_@-(u1_dddiv_@-q1_dddiv_@);r2_dddiv_@=q1_dddiv_@-r1_dddiv_@;T1L_dddiv_@=((p1_dddiv_@*r1_dddiv_@-T1H_dddiv_@)+p1_dddiv_@*r2_dddiv_@+p2_dddiv_@*r1_dddiv_@)+p2_dddiv_@*r2_dddiv_@;T1L_dddiv_@=T1L_dddiv_@+(X2L_dddiv_@*q1_dddiv_@);s1_dddiv_@=T1H_dddiv_@+T1L_dddiv_@;T1L_dddiv_@=T1L_dddiv_@-(s1_dddiv_@-T1H_dddiv_@);T1H_dddiv_@=s1_dddiv_@;s1_dddiv_@=X1H_dddiv_@-T1H_dddiv_@;r1_dddiv_@=s1_dddiv_@-X1H_dddiv_@;s2_dddiv_@=(X1H_dddiv_@-(s1_dddiv_@-r1_dddiv_@))-(T1H_dddiv_@+r1_dddiv_@);s2_dddiv_@=s2_dddiv_@-T1L_dddiv_@;s2_dddiv_@=s2_dddiv_@+X1L_dddiv_@;q2_dddiv_@=(s1_dddiv_@+s2_dddiv_@)/X2H_dddiv_@;T1H_dddiv_@=q1_dddiv_@+q2_dddiv_@;T1L_dddiv_@=q2_dddiv_@-(T1H_dddiv_@-q1_dddiv_@);c_hi=T1H_dddiv_@;c_lo=T1L_dddiv_@

#define dddiv_calc1(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_dddiv_1=a_hi;X1L_dddiv_1=a_lo;X2H_dddiv_1=b_hi;X2L_dddiv_1=b_lo;q1_dddiv_1=X1H_dddiv_1/X2H_dddiv_1;T1H_dddiv_1=X2H_dddiv_1*q1_dddiv_1;u1_dddiv_1=134217729.0d0*X2H_dddiv_1;p1_dddiv_1=u1_dddiv_1-(u1_dddiv_1-X2H_dddiv_1);p2_dddiv_1=X2H_dddiv_1-p1_dddiv_1;u1_dddiv_1=134217729.0d0*q1_dddiv_1;r1_dddiv_1=u1_dddiv_1-(u1_dddiv_1-q1_dddiv_1);r2_dddiv_1=q1_dddiv_1-r1_dddiv_1;T1L_dddiv_1=((p1_dddiv_1*r1_dddiv_1-T1H_dddiv_1)+p1_dddiv_1*r2_dddiv_1+p2_dddiv_1*r1_dddiv_1)+p2_dddiv_1*r2_dddiv_1;T1L_dddiv_1=T1L_dddiv_1+(X2L_dddiv_1*q1_dddiv_1);s1_dddiv_1=T1H_dddiv_1+T1L_dddiv_1;T1L_dddiv_1=T1L_dddiv_1-(s1_dddiv_1-T1H_dddiv_1);T1H_dddiv_1=s1_dddiv_1;s1_dddiv_1=X1H_dddiv_1-T1H_dddiv_1;r1_dddiv_1=s1_dddiv_1-X1H_dddiv_1;s2_dddiv_1=(X1H_dddiv_1-(s1_dddiv_1-r1_dddiv_1))-(T1H_dddiv_1+r1_dddiv_1);s2_dddiv_1=s2_dddiv_1-T1L_dddiv_1;s2_dddiv_1=s2_dddiv_1+X1L_dddiv_1;q2_dddiv_1=(s1_dddiv_1+s2_dddiv_1)/X2H_dddiv_1;T1H_dddiv_1=q1_dddiv_1+q2_dddiv_1;T1L_dddiv_1=q2_dddiv_1-(T1H_dddiv_1-q1_dddiv_1);c_hi=T1H_dddiv_1;c_lo=T1L_dddiv_1
#define dddiv_calc2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_dddiv_2=a_hi;X1L_dddiv_2=a_lo;X2H_dddiv_2=b_hi;X2L_dddiv_2=b_lo;q1_dddiv_2=X1H_dddiv_2/X2H_dddiv_2;T1H_dddiv_2=X2H_dddiv_2*q1_dddiv_2;u1_dddiv_2=134217729.0d0*X2H_dddiv_2;p1_dddiv_2=u1_dddiv_2-(u1_dddiv_2-X2H_dddiv_2);p2_dddiv_2=X2H_dddiv_2-p1_dddiv_2;u1_dddiv_2=134217729.0d0*q1_dddiv_2;r1_dddiv_2=u1_dddiv_2-(u1_dddiv_2-q1_dddiv_2);r2_dddiv_2=q1_dddiv_2-r1_dddiv_2;T1L_dddiv_2=((p1_dddiv_2*r1_dddiv_2-T1H_dddiv_2)+p1_dddiv_2*r2_dddiv_2+p2_dddiv_2*r1_dddiv_2)+p2_dddiv_2*r2_dddiv_2;T1L_dddiv_2=T1L_dddiv_2+(X2L_dddiv_2*q1_dddiv_2);s1_dddiv_2=T1H_dddiv_2+T1L_dddiv_2;T1L_dddiv_2=T1L_dddiv_2-(s1_dddiv_2-T1H_dddiv_2);T1H_dddiv_2=s1_dddiv_2;s1_dddiv_2=X1H_dddiv_2-T1H_dddiv_2;r1_dddiv_2=s1_dddiv_2-X1H_dddiv_2;s2_dddiv_2=(X1H_dddiv_2-(s1_dddiv_2-r1_dddiv_2))-(T1H_dddiv_2+r1_dddiv_2);s2_dddiv_2=s2_dddiv_2-T1L_dddiv_2;s2_dddiv_2=s2_dddiv_2+X1L_dddiv_2;q2_dddiv_2=(s1_dddiv_2+s2_dddiv_2)/X2H_dddiv_2;T1H_dddiv_2=q1_dddiv_2+q2_dddiv_2;T1L_dddiv_2=q2_dddiv_2-(T1H_dddiv_2-q1_dddiv_2);c_hi=T1H_dddiv_2;c_lo=T1L_dddiv_2
#define dddiv_calc3(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_dddiv_3=a_hi;X1L_dddiv_3=a_lo;X2H_dddiv_3=b_hi;X2L_dddiv_3=b_lo;q1_dddiv_3=X1H_dddiv_3/X2H_dddiv_3;T1H_dddiv_3=X2H_dddiv_3*q1_dddiv_3;u1_dddiv_3=134217729.0d0*X2H_dddiv_3;p1_dddiv_3=u1_dddiv_3-(u1_dddiv_3-X2H_dddiv_3);p2_dddiv_3=X2H_dddiv_3-p1_dddiv_3;u1_dddiv_3=134217729.0d0*q1_dddiv_3;r1_dddiv_3=u1_dddiv_3-(u1_dddiv_3-q1_dddiv_3);r2_dddiv_3=q1_dddiv_3-r1_dddiv_3;T1L_dddiv_3=((p1_dddiv_3*r1_dddiv_3-T1H_dddiv_3)+p1_dddiv_3*r2_dddiv_3+p2_dddiv_3*r1_dddiv_3)+p2_dddiv_3*r2_dddiv_3;T1L_dddiv_3=T1L_dddiv_3+(X2L_dddiv_3*q1_dddiv_3);s1_dddiv_3=T1H_dddiv_3+T1L_dddiv_3;T1L_dddiv_3=T1L_dddiv_3-(s1_dddiv_3-T1H_dddiv_3);T1H_dddiv_3=s1_dddiv_3;s1_dddiv_3=X1H_dddiv_3-T1H_dddiv_3;r1_dddiv_3=s1_dddiv_3-X1H_dddiv_3;s2_dddiv_3=(X1H_dddiv_3-(s1_dddiv_3-r1_dddiv_3))-(T1H_dddiv_3+r1_dddiv_3);s2_dddiv_3=s2_dddiv_3-T1L_dddiv_3;s2_dddiv_3=s2_dddiv_3+X1L_dddiv_3;q2_dddiv_3=(s1_dddiv_3+s2_dddiv_3)/X2H_dddiv_3;T1H_dddiv_3=q1_dddiv_3+q2_dddiv_3;T1L_dddiv_3=q2_dddiv_3-(T1H_dddiv_3-q1_dddiv_3);c_hi=T1H_dddiv_3;c_lo=T1L_dddiv_3
#define dddiv_calc4(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_dddiv_4=a_hi;X1L_dddiv_4=a_lo;X2H_dddiv_4=b_hi;X2L_dddiv_4=b_lo;q1_dddiv_4=X1H_dddiv_4/X2H_dddiv_4;T1H_dddiv_4=X2H_dddiv_4*q1_dddiv_4;u1_dddiv_4=134217729.0d0*X2H_dddiv_4;p1_dddiv_4=u1_dddiv_4-(u1_dddiv_4-X2H_dddiv_4);p2_dddiv_4=X2H_dddiv_4-p1_dddiv_4;u1_dddiv_4=134217729.0d0*q1_dddiv_4;r1_dddiv_4=u1_dddiv_4-(u1_dddiv_4-q1_dddiv_4);r2_dddiv_4=q1_dddiv_4-r1_dddiv_4;T1L_dddiv_4=((p1_dddiv_4*r1_dddiv_4-T1H_dddiv_4)+p1_dddiv_4*r2_dddiv_4+p2_dddiv_4*r1_dddiv_4)+p2_dddiv_4*r2_dddiv_4;T1L_dddiv_4=T1L_dddiv_4+(X2L_dddiv_4*q1_dddiv_4);s1_dddiv_4=T1H_dddiv_4+T1L_dddiv_4;T1L_dddiv_4=T1L_dddiv_4-(s1_dddiv_4-T1H_dddiv_4);T1H_dddiv_4=s1_dddiv_4;s1_dddiv_4=X1H_dddiv_4-T1H_dddiv_4;r1_dddiv_4=s1_dddiv_4-X1H_dddiv_4;s2_dddiv_4=(X1H_dddiv_4-(s1_dddiv_4-r1_dddiv_4))-(T1H_dddiv_4+r1_dddiv_4);s2_dddiv_4=s2_dddiv_4-T1L_dddiv_4;s2_dddiv_4=s2_dddiv_4+X1L_dddiv_4;q2_dddiv_4=(s1_dddiv_4+s2_dddiv_4)/X2H_dddiv_4;T1H_dddiv_4=q1_dddiv_4+q2_dddiv_4;T1L_dddiv_4=q2_dddiv_4-(T1H_dddiv_4-q1_dddiv_4);c_hi=T1H_dddiv_4;c_lo=T1L_dddiv_4
#define dddiv_calc5(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_dddiv_5=a_hi;X1L_dddiv_5=a_lo;X2H_dddiv_5=b_hi;X2L_dddiv_5=b_lo;q1_dddiv_5=X1H_dddiv_5/X2H_dddiv_5;T1H_dddiv_5=X2H_dddiv_5*q1_dddiv_5;u1_dddiv_5=134217729.0d0*X2H_dddiv_5;p1_dddiv_5=u1_dddiv_5-(u1_dddiv_5-X2H_dddiv_5);p2_dddiv_5=X2H_dddiv_5-p1_dddiv_5;u1_dddiv_5=134217729.0d0*q1_dddiv_5;r1_dddiv_5=u1_dddiv_5-(u1_dddiv_5-q1_dddiv_5);r2_dddiv_5=q1_dddiv_5-r1_dddiv_5;T1L_dddiv_5=((p1_dddiv_5*r1_dddiv_5-T1H_dddiv_5)+p1_dddiv_5*r2_dddiv_5+p2_dddiv_5*r1_dddiv_5)+p2_dddiv_5*r2_dddiv_5;T1L_dddiv_5=T1L_dddiv_5+(X2L_dddiv_5*q1_dddiv_5);s1_dddiv_5=T1H_dddiv_5+T1L_dddiv_5;T1L_dddiv_5=T1L_dddiv_5-(s1_dddiv_5-T1H_dddiv_5);T1H_dddiv_5=s1_dddiv_5;s1_dddiv_5=X1H_dddiv_5-T1H_dddiv_5;r1_dddiv_5=s1_dddiv_5-X1H_dddiv_5;s2_dddiv_5=(X1H_dddiv_5-(s1_dddiv_5-r1_dddiv_5))-(T1H_dddiv_5+r1_dddiv_5);s2_dddiv_5=s2_dddiv_5-T1L_dddiv_5;s2_dddiv_5=s2_dddiv_5+X1L_dddiv_5;q2_dddiv_5=(s1_dddiv_5+s2_dddiv_5)/X2H_dddiv_5;T1H_dddiv_5=q1_dddiv_5+q2_dddiv_5;T1L_dddiv_5=q2_dddiv_5-(T1H_dddiv_5-q1_dddiv_5);c_hi=T1H_dddiv_5;c_lo=T1L_dddiv_5
