  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! 
  !-----------------------------------------------------------------------
  SUBROUTINE fermiwindow 
  !-----------------------------------------------------------------------
  !
  !  find the band indices of the first
  !  and last state falling within the window e_fermi+-efermithickness
  ! 
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE el_phon,       ONLY : etf, ibndmin, ibndmax, nksf
  USE epwcom,        ONLY : fsthick, fsthickprime, efprime, nbndsub, egap
  USE pwcom,         ONLY : ef
  USE control_flags, ONLY : iverbosity
#ifdef __PARA
  USE mp,            ONLY : mp_max, mp_min
  USE mp_global,     ONLY : inter_pool_comm
#endif
  implicit none
  integer :: ik, ibnd
  real(kind=DP) :: tmp, ebnd, ebndmin, ebndmax,emin_above_ef, emax_below_ef
  integer :: ibnd_max_below_ef, ibnd_min_above_ef, is_metal
  !
  ibndmin = 100000
  ibndmax = 0
  ebndmin =  1.d8
  ebndmax = -1.d8

  emin_above_ef= 1.d8
  emax_below_ef= -1.d8
  ibnd_max_below_ef=0
  ibnd_min_above_ef=10000
  is_metal=0
  !
  DO ik = 1, nksf
    DO ibnd = 1, nbndsub
      ebnd = etf (ibnd, ik)
      !
      IF  ( abs(ebnd - efprime) .lt. fsthickprime ) THEN
      ! the inclusion of bands should consider all choice of efs
        ibndmin = min(ibnd,ibndmin)
        ibndmax = max(ibnd,ibndmax)
        ebndmin = min(ebnd,ebndmin)
        ebndmax = max(ebnd,ebndmax)
      ENDIF
      !Finding band gap
      IF (ebnd .le. ef .and. ebnd .gt. emax_below_ef) THEN
         emax_below_ef = ebnd
         ibnd_max_below_ef = ibnd
      ENDIF
      IF (ebnd .ge. ef .and. ebnd .lt. emin_above_ef) THEN
         emin_above_ef = ebnd
         ibnd_min_above_ef = ibnd
      ENDIF
      IF (ibnd_min_above_ef .le. ibnd_max_below_ef) THEN
         is_metal=1
      ENDIF
      !
    ENDDO
  ENDDO
  !
#ifdef __PARA
  tmp = float (ibndmin)
  CALL mp_min(tmp,inter_pool_comm)
  ibndmin = nint (tmp)
 CALL mp_min(ebndmin,inter_pool_comm)
  !
  tmp = float (ibndmax)
  CALL mp_max(tmp, inter_pool_comm)
  ibndmax = nint (tmp)
  CALL mp_max(ebndmax,inter_pool_comm)
  !
  tmp = emin_above_ef
  CALL mp_min(tmp, inter_pool_comm)
  emin_above_ef = tmp
  !
  tmp = emax_below_ef
  CALL mp_max(tmp, inter_pool_comm)
  emax_below_ef = tmp
  !
  tmp = float (is_metal)
  CALL mp_max(tmp, inter_pool_comm)
  is_metal = nint(tmp)
#endif
  !
  IF (is_metal .eq. 1) THEN
       egap = 0
       WRITE(stdout,'(/14x,a)') 'This is a metal, egap = 0'
  ELSE
       egap = emin_above_ef - emax_below_ef
       WRITE(stdout,'(/14x,a,f9.3,a)') 'Energy gap = ',egap*13.60569,'eV'
  ENDIF 
 ! IF (iverbosity.eq.1) THEN
    WRITE(stdout,'(/14x,a,i5,2x,a,f9.3)') 'ibndmin = ', ibndmin, 'ebndmin = ', ebndmin
    WRITE(stdout,'(14x,a,i5,2x,a,f9.3/)') 'ibndmax = ', ibndmax, 'ebndmax = ', ebndmax
 ! ENDIF
  !
  END SUBROUTINE fermiwindow

