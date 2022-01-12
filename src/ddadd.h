! -------------
! add
! real*8 T1H,T1L
! real*8 X1H,X1L
! real*8 X2H,X2L
! real*8 p1,p2
! real*8 s1

! #define ddadd_var@ T1H_ddadd_@,T1L_ddadd_@,X1H_ddadd_@,X1L_ddadd_@,X2H_ddadd_@,X2L_ddadd_@,p1_ddadd_@,p2_ddadd_@,s1
#define ddadd_var1 T1H_ddadd_1, T1L_ddadd_1, X1H_ddadd_1, X1L_ddadd_1, X2H_ddadd_1, X2L_ddadd_1, p1_ddadd_1, p2_ddadd_1, s1_ddadd_1
#define ddadd_var2 T1H_ddadd_2, T1L_ddadd_2, X1H_ddadd_2, X1L_ddadd_2, X2H_ddadd_2, X2L_ddadd_2, p1_ddadd_2, p2_ddadd_2, s1_ddadd_2
#define ddadd_var3 T1H_ddadd_3, T1L_ddadd_3, X1H_ddadd_3, X1L_ddadd_3, X2H_ddadd_3, X2L_ddadd_3, p1_ddadd_3, p2_ddadd_3, s1_ddadd_3
#define ddadd_var4 T1H_ddadd_4, T1L_ddadd_4, X1H_ddadd_4, X1L_ddadd_4, X2H_ddadd_4, X2L_ddadd_4, p1_ddadd_4, p2_ddadd_4, s1_ddadd_4
#define ddadd_var5 T1H_ddadd_5, T1L_ddadd_5, X1H_ddadd_5, X1L_ddadd_5, X2H_ddadd_5, X2L_ddadd_5, p1_ddadd_5, p2_ddadd_5, s1_ddadd_5

! -------------    
! X1H_ddadd_@  = a_hi
! X1L_ddadd_@  = a_lo
! X2H_ddadd_@  = b_hi
! X2L_ddadd_@  = b_lo
! p1_ddadd_@   = X1H_ddadd_@ + X2H_ddadd_@
! p2_ddadd_@   = p1_ddadd_@ - X1H_ddadd_@
! T1L_ddadd_@  = (X1H_ddadd_@ - (p1_ddadd_@ - p2_ddadd_@)) + (X2H_ddadd_@ - p2_ddadd_@)
! T1H_ddadd_@  = p1_ddadd_@
! T1L_ddadd_@  = T1L_ddadd_@ + X1L_ddadd_@ + X2L_ddadd_@
! s1_ddadd_@   = T1H_ddadd_@ + T1L_ddadd_@
! T1L_ddadd_@  = T1L_ddadd_@ - (s1_ddadd_@ - T1H_ddadd_@)
! T1H_ddadd_@  = s1_ddadd_@
! c_hi = T1H_ddadd_@
! c_lo = T1L_ddadd_@
! #define ddadd_calc(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddadd_@=a_hi;X1L_ddadd_@=a_lo;X2H_ddadd_@=b_hi;X2L_ddadd_@=b_lo;p1_ddadd_@=X1H_ddadd_@+X2H_ddadd_@;p2_ddadd_@=p1_ddadd_@-X1H_ddadd_@;T1L_ddadd_@=(X1H_ddadd_@-(p1_ddadd_@-p2_ddadd_@))+(X2H_ddadd_@-p2_ddadd_@);T1H_ddadd_@=p1_ddadd_@;T1L_ddadd_@=T1L_ddadd_@+X1L_ddadd_@+X2L_ddadd_@;s1_ddadd_@=T1H_ddadd_@+T1L_ddadd_@;T1L_ddadd_@=T1L_ddadd_@-(s1_ddadd_@-T1H_ddadd_@);T1H_ddadd_@=s1_ddadd_@;c_hi=T1H_ddadd_@;c_lo=T1L_ddadd_@

#define ddadd_calc1(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddadd_1=a_hi;X1L_ddadd_1=a_lo;X2H_ddadd_1=b_hi;X2L_ddadd_1=b_lo;p1_ddadd_1=X1H_ddadd_1+X2H_ddadd_1;p2_ddadd_1=p1_ddadd_1-X1H_ddadd_1;T1L_ddadd_1=(X1H_ddadd_1-(p1_ddadd_1-p2_ddadd_1))+(X2H_ddadd_1-p2_ddadd_1);T1H_ddadd_1=p1_ddadd_1;T1L_ddadd_1=T1L_ddadd_1+X1L_ddadd_1+X2L_ddadd_1;s1_ddadd_1=T1H_ddadd_1+T1L_ddadd_1;T1L_ddadd_1=T1L_ddadd_1-(s1_ddadd_1-T1H_ddadd_1);T1H_ddadd_1=s1_ddadd_1;c_hi=T1H_ddadd_1;c_lo=T1L_ddadd_1

#define ddadd_calc2(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddadd_2=a_hi;X1L_ddadd_2=a_lo;X2H_ddadd_2=b_hi;X2L_ddadd_2=b_lo;p1_ddadd_2=X1H_ddadd_2+X2H_ddadd_2;p2_ddadd_2=p1_ddadd_2-X1H_ddadd_2;T1L_ddadd_2=(X1H_ddadd_2-(p1_ddadd_2-p2_ddadd_2))+(X2H_ddadd_2-p2_ddadd_2);T1H_ddadd_2=p1_ddadd_2;T1L_ddadd_2=T1L_ddadd_2+X1L_ddadd_2+X2L_ddadd_2;s1_ddadd_2=T1H_ddadd_2+T1L_ddadd_2;T1L_ddadd_2=T1L_ddadd_2-(s1_ddadd_2-T1H_ddadd_2);T1H_ddadd_2=s1_ddadd_2;c_hi=T1H_ddadd_2;c_lo=T1L_ddadd_2

#define ddadd_calc3(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddadd_3=a_hi;X1L_ddadd_3=a_lo;X2H_ddadd_3=b_hi;X2L_ddadd_3=b_lo;p1_ddadd_3=X1H_ddadd_3+X2H_ddadd_3;p2_ddadd_3=p1_ddadd_3-X1H_ddadd_3;T1L_ddadd_3=(X1H_ddadd_3-(p1_ddadd_3-p2_ddadd_3))+(X2H_ddadd_3-p2_ddadd_3);T1H_ddadd_3=p1_ddadd_3;T1L_ddadd_3=T1L_ddadd_3+X1L_ddadd_3+X2L_ddadd_3;s1_ddadd_3=T1H_ddadd_3+T1L_ddadd_3;T1L_ddadd_3=T1L_ddadd_3-(s1_ddadd_3-T1H_ddadd_3);T1H_ddadd_3=s1_ddadd_3;c_hi=T1H_ddadd_3;c_lo=T1L_ddadd_3

#define ddadd_calc4(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddadd_4=a_hi;X1L_ddadd_4=a_lo;X2H_ddadd_4=b_hi;X2L_ddadd_4=b_lo;p1_ddadd_4=X1H_ddadd_4+X2H_ddadd_4;p2_ddadd_4=p1_ddadd_4-X1H_ddadd_4;T1L_ddadd_4=(X1H_ddadd_4-(p1_ddadd_4-p2_ddadd_4))+(X2H_ddadd_4-p2_ddadd_4);T1H_ddadd_4=p1_ddadd_4;T1L_ddadd_4=T1L_ddadd_4+X1L_ddadd_4+X2L_ddadd_4;s1_ddadd_4=T1H_ddadd_4+T1L_ddadd_4;T1L_ddadd_4=T1L_ddadd_4-(s1_ddadd_4-T1H_ddadd_4);T1H_ddadd_4=s1_ddadd_4;c_hi=T1H_ddadd_4;c_lo=T1L_ddadd_4

#define ddadd_calc5(a_hi,a_lo,b_hi,b_lo,c_hi,c_lo) X1H_ddadd_5=a_hi;X1L_ddadd_5=a_lo;X2H_ddadd_5=b_hi;X2L_ddadd_5=b_lo;p1_ddadd_5=X1H_ddadd_5+X2H_ddadd_5;p2_ddadd_5=p1_ddadd_5-X1H_ddadd_5;T1L_ddadd_5=(X1H_ddadd_5-(p1_ddadd_5-p2_ddadd_5))+(X2H_ddadd_5-p2_ddadd_5);T1H_ddadd_5=p1_ddadd_5;T1L_ddadd_5=T1L_ddadd_5+X1L_ddadd_5+X2L_ddadd_5;s1_ddadd_5=T1H_ddadd_5+T1L_ddadd_5;T1L_ddadd_5=T1L_ddadd_5-(s1_ddadd_5-T1H_ddadd_5);T1H_ddadd_5=s1_ddadd_5;c_hi=T1H_ddadd_5;c_lo=T1L_ddadd_5
