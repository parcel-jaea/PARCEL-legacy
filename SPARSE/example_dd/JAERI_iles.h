!
! header file for iterative linear equation solvers
!

! control parameters
! max tolerable value of the relative residual norm
      real*8  Res_max
! order of Neumann expansion for additive Schwartz preconditioner
      integer  N_Neumann

      integer AS_type
      integer ILU_method 
      common /jaeri_param/ Res_max, N_Neumann, AS_type, ILU_method
      save /jaeri_param/

! debug out
      integer  dbgout
      common /jaeri_dbg/ dbgout

! timing
      double precision  tw(6), tc(6)
      common /jaeri_time/  tw, tc

! standard error output
      integer  stderr
      parameter (stderr = 0)

! return code
      integer  RETURN_SUCCESS, RETURN_NOT_CONVERGENT, RETURN_DIVERGENT
      parameter( RETURN_SUCCESS        = 0 )
      parameter( RETURN_NOT_CONVERGENT = 1 )
      parameter( RETURN_DIVERGENT      = 2)

! factorization metod in local(additive schwarz) block
!
! LF_ILUD: ILUD(D-ILU)
! LF_MILUD: MILUD
!
      integer  LF_ILUD, LF_MILUD
      parameter (LF_ILUD  = 1, LF_MILUD = 2)

! type of interpolation for the overlaps of preconditioned vectors
      integer  AS_NO_INTERPOLATE, AS_INTERPOLATE
      parameter (AS_NO_INTERPOLATE = 0, AS_INTERPOLATE = 1)

