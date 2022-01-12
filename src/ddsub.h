  ! -------------
  ! real*8 T1H,T1L
  ! real*8 X1H,X1L
  ! real*8 X2H,X2L
  ! real*8 p1,p2
  ! real*8 s1
  ! #define ddsub_var@ T1H_ddsub_@,T1L_ddsub_@,X1H_ddsub_@,X1L_ddsub_@,X2H_ddsub_@,X2L_ddsub_@,p1_ddsub_@,p2_ddsub_@,s1_ddsub_@
#define ddsub_var1 T1H_ddsub_1,T1L_ddsub_1,X1H_ddsub_1,X1L_ddsub_1,X2H_ddsub_1,X2L_ddsub_1,p1_ddsub_1,p2_ddsub_1,s1_ddsub_1
#define ddsub_var2 T1H_ddsub_2,T1L_ddsub_2,X1H_ddsub_2,X1L_ddsub_2,X2H_ddsub_2,X2L_ddsub_2,p1_ddsub_2,p2_ddsub_2,s1_ddsub_2
#define ddsub_var3 T1H_ddsub_3,T1L_ddsub_3,X1H_ddsub_3,X1L_ddsub_3,X2H_ddsub_3,X2L_ddsub_3,p1_ddsub_3,p2_ddsub_3,s1_ddsub_3
#define ddsub_var4 T1H_ddsub_4,T1L_ddsub_4,X1H_ddsub_4,X1L_ddsub_4,X2H_ddsub_4,X2L_ddsub_4,p1_ddsub_4,p2_ddsub_4,s1_ddsub_4
#define ddsub_var5 T1H_ddsub_5,T1L_ddsub_5,X1H_ddsub_5,X1L_ddsub_5,X2H_ddsub_5,X2L_ddsub_5,p1_ddsub_5,p2_ddsub_5,s1_ddsub_5

  ! ------------
    ! X1H_ddsub_@  = a_hi
    ! X1L_ddsub_@  = a_lo
    ! X2H_ddsub_@  = b_hi
    ! X2L_ddsub_@  = b_lo
    ! p1_ddsub_@   = X1H_ddsub_@ - X2H_ddsub_@
    ! p2_ddsub_@   = p1_ddsub_@ - X1H_ddsub_@
    ! T1L_ddsub_@  = (X1H_ddsub_@ - (p1_ddsub_@ - p2_ddsub_@)) - (X2H_ddsub_@ + p2_ddsub_@)
    ! T1H_ddsub_@  = p1_ddsub_@
    ! T1L_ddsub_@  = T1L_ddsub_@ + X1L_ddsub_@ - X2L_ddsub_@
    ! s1_ddsub_@   = T1H_ddsub_@ + T1L_ddsub_@
    ! T1L_ddsub_@  = T1L_ddsub_@ - (s1_ddsub_@ - T1H_ddsub_@)
    ! T1H_ddsub_@  = s1_ddsub_@
    ! c_hi = T1H_ddsub_@
    ! c_lo = T1L_ddsub_@

! #define ddsub_calc@(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddsub_@=a_hi;X1L_ddsub_@=a_lo;X2H_ddsub_@=b_hi;X2L_ddsub_@=b_lo;p1_ddsub_@=X1H_ddsub_@-X2H_ddsub_@;p2_ddsub_@=p1_ddsub_@-X1H_ddsub_@;T1L_ddsub_@=(X1H_ddsub_@-(p1_ddsub_@-p2_ddsub_@))-(X2H_ddsub_@+p2_ddsub_@);T1H_ddsub_@=p1_ddsub_@;T1L_ddsub_@=T1L_ddsub_@+X1L_ddsub_@-X2L_ddsub_@;s1_ddsub_@=T1H_ddsub_@+T1L_ddsub_@;T1L_ddsub_@=T1L_ddsub_@-(s1_ddsub_@-T1H_ddsub_@);T1H_ddsub_@=s1_ddsub_@;c_hi=T1H_ddsub_@;c_lo=T1L_ddsub_@

#define ddsub_calc1(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddsub_1=a_hi;X1L_ddsub_1=a_lo;X2H_ddsub_1=b_hi;X2L_ddsub_1=b_lo;p1_ddsub_1=X1H_ddsub_1-X2H_ddsub_1;p2_ddsub_1=p1_ddsub_1-X1H_ddsub_1;T1L_ddsub_1=(X1H_ddsub_1-(p1_ddsub_1-p2_ddsub_1))-(X2H_ddsub_1+p2_ddsub_1);T1H_ddsub_1=p1_ddsub_1;T1L_ddsub_1=T1L_ddsub_1+X1L_ddsub_1-X2L_ddsub_1;s1_ddsub_1=T1H_ddsub_1+T1L_ddsub_1;T1L_ddsub_1=T1L_ddsub_1-(s1_ddsub_1-T1H_ddsub_1);T1H_ddsub_1=s1_ddsub_1;c_hi=T1H_ddsub_1;c_lo=T1L_ddsub_1
#define ddsub_calc2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddsub_2=a_hi;X1L_ddsub_2=a_lo;X2H_ddsub_2=b_hi;X2L_ddsub_2=b_lo;p1_ddsub_2=X1H_ddsub_2-X2H_ddsub_2;p2_ddsub_2=p1_ddsub_2-X1H_ddsub_2;T1L_ddsub_2=(X1H_ddsub_2-(p1_ddsub_2-p2_ddsub_2))-(X2H_ddsub_2+p2_ddsub_2);T1H_ddsub_2=p1_ddsub_2;T1L_ddsub_2=T1L_ddsub_2+X1L_ddsub_2-X2L_ddsub_2;s1_ddsub_2=T1H_ddsub_2+T1L_ddsub_2;T1L_ddsub_2=T1L_ddsub_2-(s1_ddsub_2-T1H_ddsub_2);T1H_ddsub_2=s1_ddsub_2;c_hi=T1H_ddsub_2;c_lo=T1L_ddsub_2
#define ddsub_calc3(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddsub_3=a_hi;X1L_ddsub_3=a_lo;X2H_ddsub_3=b_hi;X2L_ddsub_3=b_lo;p1_ddsub_3=X1H_ddsub_3-X2H_ddsub_3;p2_ddsub_3=p1_ddsub_3-X1H_ddsub_3;T1L_ddsub_3=(X1H_ddsub_3-(p1_ddsub_3-p2_ddsub_3))-(X2H_ddsub_3+p2_ddsub_3);T1H_ddsub_3=p1_ddsub_3;T1L_ddsub_3=T1L_ddsub_3+X1L_ddsub_3-X2L_ddsub_3;s1_ddsub_3=T1H_ddsub_3+T1L_ddsub_3;T1L_ddsub_3=T1L_ddsub_3-(s1_ddsub_3-T1H_ddsub_3);T1H_ddsub_3=s1_ddsub_3;c_hi=T1H_ddsub_3;c_lo=T1L_ddsub_3
#define ddsub_calc4(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddsub_4=a_hi;X1L_ddsub_4=a_lo;X2H_ddsub_4=b_hi;X2L_ddsub_4=b_lo;p1_ddsub_4=X1H_ddsub_4-X2H_ddsub_4;p2_ddsub_4=p1_ddsub_4-X1H_ddsub_4;T1L_ddsub_4=(X1H_ddsub_4-(p1_ddsub_4-p2_ddsub_4))-(X2H_ddsub_4+p2_ddsub_4);T1H_ddsub_4=p1_ddsub_4;T1L_ddsub_4=T1L_ddsub_4+X1L_ddsub_4-X2L_ddsub_4;s1_ddsub_4=T1H_ddsub_4+T1L_ddsub_4;T1L_ddsub_4=T1L_ddsub_4-(s1_ddsub_4-T1H_ddsub_4);T1H_ddsub_4=s1_ddsub_4;c_hi=T1H_ddsub_4;c_lo=T1L_ddsub_4
#define ddsub_calc5(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddsub_5=a_hi;X1L_ddsub_5=a_lo;X2H_ddsub_5=b_hi;X2L_ddsub_5=b_lo;p1_ddsub_5=X1H_ddsub_5-X2H_ddsub_5;p2_ddsub_5=p1_ddsub_5-X1H_ddsub_5;T1L_ddsub_5=(X1H_ddsub_5-(p1_ddsub_5-p2_ddsub_5))-(X2H_ddsub_5+p2_ddsub_5);T1H_ddsub_5=p1_ddsub_5;T1L_ddsub_5=T1L_ddsub_5+X1L_ddsub_5-X2L_ddsub_5;s1_ddsub_5=T1H_ddsub_5+T1L_ddsub_5;T1L_ddsub_5=T1L_ddsub_5-(s1_ddsub_5-T1H_ddsub_5);T1H_ddsub_5=s1_ddsub_5;c_hi=T1H_ddsub_5;c_lo=T1L_ddsub_5
