/*******************************************************************
 *
 *    Elapsed CPU time monitor
 *
 *    Programmed by M.Sasaki (the Japan Research Institute Ltd.)
 *
 *    Used with mutually exclusive cpp options;
 *       -DPOSIX_C, -DSUN_C  and -DANSI_C.
 *
 ******************************************************************/
/******************************************************************
 *  Fortran callable C subroutine :
 *
 *   UNIX venders have their own convention or rule for external
 *  names of FORTRAN coded subprograms or commons.
 *
 *   We prepare equivalent sources for C coded routines
 * 
 *     1. Attach underscore to FORTRAN name. ( SUN & Other BSD UNIX)
 *     2. Same as FORTRAN name (lowercase).   (HP-UX,IBM-AIX etc.)
 *     3. Uppercase string of FORTRAN name.   ( Ncube etc. )
 *     4. "_" + uppercase string of FORTRAN name.   ( HI/OSF )
 *
 *******************************************************************/

#include <stdio.h>

#ifdef SUN_C

#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>

#elif POSIX_C 

#include <unistd.h>
#include <sys/times.h>

#elif ANSI_C 

#include <time.h>

#endif



#ifdef POSIX_C

/*******************************************************************
 *
 *    For environments having POSIX library routine "times".
 *    ( the function "times" is defined in POSIX but not in ANSI-C)
 *
 ******************************************************************/

static int init =0 ;

void cputime_( ttt )
   double * ttt ;
{
   static struct tms  tt ;
   static double resolution ;

   if( init == 0 )
   {
       resolution = 1.0/(double)sysconf( _SC_CLK_TCK ) ;
       times( &tt ) ;
       init = 1 ;
   }
   times( &tt ) ;
   /**return (tt.tms_utime + tt.tms_stime)*resolution **/;
   *ttt =  (double)((tt.tms_utime + tt.tms_stime)*resolution) ;
}

void cputime ( ttt )
   double * ttt ;
{
   static struct tms  tt ;
   static double resolution ;

   if( init == 0 )
   {
       resolution = 1.0/(double)sysconf( _SC_CLK_TCK ) ;
       times( &tt ) ;
       init = 1 ;
   }
   times( &tt ) ;
   *ttt =  (double)((tt.tms_utime + tt.tms_stime)*resolution) ;
}
void CPUTIME( ttt )
   double * ttt ;
{
   static struct tms  tt ;
   static double resolution ;

   if( init == 0 )
   {
       resolution = 1.0/(double)sysconf( _SC_CLK_TCK ) ;
       times( &tt ) ;
       init = 1 ;
   }
   times( &tt ) ;
   *ttt =  (double)((tt.tms_utime + tt.tms_stime)*resolution) ;
}
void _CPUTIME( ttt )
   double * ttt ;
{
   static struct tms  tt ;
   static double resolution ;

   if( init == 0 )
   {
       resolution = 1.0/(double)sysconf( _SC_CLK_TCK ) ;
       times( &tt ) ;
       init = 1 ;
   }
   times( &tt ) ;
   *ttt =  (double)((tt.tms_utime + tt.tms_stime)*resolution) ;
}

#endif /** end of POSIX_C version **/

#ifdef SUN_C

/*******************************************************************
 *
 *    Elapsed CPU time monitor on old SUN-OS's whose clock tic
 *   may be 1/60 sec.
 *
 *    Programmed by M.Sasaki (the Japan Research Institute Co. Ltd.)
 *    
 *    1st version:  20 Jan 1993
 *
 ******************************************************************/


static int init =0 ;

void cputime_( ttt )
   double * ttt ;
{
   static struct tms  tt ;
   static double resolution ;

   if( init == 0 )
   {
#ifdef _SC_CLK_TCK
       resolution = 1.0/(double)sysconf( _SC_CLK_TCK ) ;
#else
       resolution = 1.0/(double)60.0 ;
#endif
       times( &tt ) ;
       init = 1 ;
   }

   times( &tt ) ;
   *ttt =  (double)((tt.tms_utime + tt.tms_stime)*resolution) ;
}

#endif /** end of SUN_C version **/

#ifdef ANSI_C

/*******************************************************************
 *  Accumulated "CPU" (?) time using ANSI C function clock().
 *
 *  For systems whose clock tic is microsecond and clock_t is
 *  32 bit integer:
 *
 *  Calling intervals longer than 36 minutes may give inaccurate
 *  timings, and it is impossible to correct timings for calling
 *  intervals more than 72 minutes :-).
 *******************************************************************/

static int init =0 ;
void MY_CPUTIME() ;

void cputime( ttt )
   double * ttt ;
{
   MY_CPUTIME(ttt) ;
}
void cputime_( ttt )
   double * ttt ;
{
   MY_CPUTIME(ttt) ;
}
void CPUTIME( ttt )
   double * ttt ;
{
   MY_CPUTIME(ttt) ;
}
void _CPUTIME( ttt )
   double * ttt ;
{
   MY_CPUTIME(ttt) ;
}

void MY_CPUTIME( ttt )
   double * ttt ;
{
   static clock_t tt, t1;
   static double resolution, tttd ;
   static double accum_time = 0.0 ;

#define MIN_36 35.791394133333   /*** 2^31 / 1000,000 / 60 ***/

   if( init == 0 )
   {
       tt   = clock() ;
       accum_time = (double)tt / CLOCKS_PER_SEC ;

       if ( CLOCKS_PER_SEC == 1000000 )
       {
          init = 1 ;  /** clock tic is micro second  **/
       }
       else
          init = 2 ;

       resolution = 1.0/(double)CLOCKS_PER_SEC ;
   }

   t1 = clock() ;
   accum_time += (double)(t1 - tt)*resolution ;

  /******************************************************
   *  Trying to correct timings for calling intevals
   *  longer than 36 minutes.
   ******************************************************/
   if( init == 1 && t1 < tt ) accum_time += MIN_36 ;

   tt   = t1 ;
   *ttt = (double)accum_time ;
}
#endif /** end of ANSI_C version **/

