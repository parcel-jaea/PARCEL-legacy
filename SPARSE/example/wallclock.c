
/*******************************************************************
 *
 *    Wall clock time since midnight January 1, 1970. ;
 *
 *    This clock has a nominal precision of microsecond but
 *  the actual precision depends on settings of users system.
 *
 *    Programmed by M.Sasaki (the Japan Research Institute Co. Ltd.)
 *    
 *    1st version:  6 Oct 1993
 *    19 Jul 1996 : use time() for non-POSIX and AINS-C compilers.
 *
 *   * works on HP-UX, SunOS, DEC/OSF, Sysytem_V etc. 
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

#ifdef ANSI_C

 /*** Use time function for ANSI-C (non POSIX) mode ***/

#include <time.h>
time_t ttt0 ;

double wall_clock()
{
    return (double)time( &ttt0 ) ;
}

#else

 /*** POSIX:  use gettimeofday **/

#include <sys/types.h>
#include <sys/time.h>

double wall_clock()
{
     struct timeval  tv;
     struct timezone  tz;
#ifdef SX        /*** Super UX lacks yz argument of gettimeofday **/
     gettimeofday( &tv ) ;
#else
     gettimeofday( &tv, &tz ) ;
#endif
     return (double)tv.tv_sec + (double)(tv.tv_usec)/1000000 ;
}
#endif

void wallclock( ttt )
   double *ttt ;
{
    *ttt = wall_clock();
}
void wallclock_( ttt )
   double *ttt ;
{
    *ttt = wall_clock();
}
void WALLCLOCK( ttt )
   double *ttt ;
{
    *ttt = wall_clock();
}
void _WALLCLOCK( ttt )
   double *ttt ;
{
    *ttt = wall_clock();
}
