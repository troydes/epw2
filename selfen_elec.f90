  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !   
  
  
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_elec_fly ( iq )
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the electron self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the electron linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  USE phcom, ONLY : lgamma, nmodes
  USE epwcom, ONLY : nbndsub, lrepmatf, iunepmatf, &
      fsthick, eptemp, ngaussw, degaussw, iuetf,   &
      nbndskip, ecutse, parallel_k, ndeg,  &
      parallel_q, epf_mem, etf_mem, eig_read, eps_acustic, &
      nefsweep, isefsweep, neptemp, efsweep, &
      perphon, constepm, groupq, nqgroup
  USE pwcom, ONLY : nelec, ef, isk
  USE el_phon, ONLY : etf, ibndmin, ibndmax, nksf, etfq, &
      epf17, wkf, nksqf, nxqf, wf, wqf, xkf, nkstotf, &
      sigmai_all_q, ekk_all, xkf_all, sigmai_all_gpq_abs, sigmai_all_gpq_emi, gpid
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_global, ONLY : me_pool, inter_pool_comm, my_pool_id
#endif
  implicit none
  !
  REAL(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0)
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount, ideg, iefs, ieptemp
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  REAL(kind=DP) :: g2, ekk, ekq, wq, ef0(ndeg), wgq, wgkq,  &
       weight, wgauss, dosef(ndeg), dos_ef, sigmar(nbndsub, nksf, ndeg, neptemp, nefsweep), &
       sigmai(nbndsub, nksf, ndeg, neptemp, nefsweep), zi(nbndsub, nksf, ndeg, neptemp, nefsweep), &
       eptemp0, degausswi, ef1, sigmai_mode_abs(nbndsub, nksf, ndeg, neptemp, nefsweep, nmodes), &
       sigmai_mode_emi(nbndsub, nksf, ndeg, neptemp, nefsweep, nmodes), weight_abs, weight_emi
  logical :: already_skipped
  REAL(kind=DP), external :: efermig
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf
  REAL(kind=DP), allocatable :: etf_all(:,:), &
                                sigmar_all (:,:,:,:,:), sigmai_all (:,:,:,:,:), zi_all (:,:,:,:,:), &
                                sigmai_all_mode_abs (:,:,:,:,:,:), sigmai_all_mode_emi (:,:,:,:,:,:)
  !
 
  !
  ! loop over temperatures can be introduced
  !
  eptemp0 = eptemp(1)
  !
!   IF ( fsthick .lt. 1.d3 ) &
!     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
!       'Fermi Surface thickness = ', fsthick, ' Ry'
!   WRITE(stdout, '(/5x,a,e18.9,a)' ) &
!     'Golden Rule strictly enforced with T = ',eptemp0, ' Ry'
  !
  ! here we take into account that we may skip bands when we wannierize
  ! (spin-unpolarized)
  already_skipped = .false.
  IF ( nbndskip .gt. 0 ) THEN
     IF ( .not. already_skipped ) THEN 
        nelec = nelec - two * nbndskip
        already_skipped = .true.
        WRITE(stdout, '(/5x,"Skipping the first ",i4," bands:")' ) nbndskip
        WRITE(stdout, '(/5x,"The Fermi level will be determined with ",f9.5," electrons")' ) nelec
     ENDIF
  ENDIF

  !
  ! Fermi level and corresponding DOS
  !
  if ( isefsweep .eq. 0) then !! when sweeping through efs, the dosef(ideg) below will be wrong. So now phonon_self cannot be done
          do ideg = 1, ndeg
                ef0(ideg) = efermig(etf,nbndsub,nksf,nelec,wkf,degaussw(ideg),ngaussw,0,isk)
                  !
                  !   if 'fine' Fermi level differs by more than 250 meV, there is probably
                  !   something wrong with the wannier functions
                  IF ( abs(ef0(ideg) - ef) * ryd2eV .gt. 0.5 .and. (.not.eig_read) ) &
                       CALL errore ('selfen_elec', 'Something wrong with Ef, check MLWFs', 1)
                  !
                dosef(ideg) = dos_ef (ngaussw, degaussw(ideg), ef0(ideg), etf, wkf, nksf, nbndsub)
                !   N(Ef) in the equation for lambda is the DOS per spin
                dosef(ideg) = dosef(ideg) / two
          enddo
        endif
  
  
  !
!   WRITE (6, 100) degaussw, ngaussw
!   WRITE (6, 101) dosef, ef0 * ryd2ev
!   WRITE (6, 101) dosef, ef  * ryd2ev
!   WRITE (6,'(a)') ' '
  !
  sigmar = zero
  sigmai = zero
  zi = zero
  sigmai_mode_abs = zero
  sigmai_mode_emi = zero
  !
  ! loop over all k points of the fine mesh
  !
  fermicount = 0 
  DO ik = 1, nksqf
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     !
     ! No loop over the q-points
     !
!!     DO iq = 1, nxqf
        !
        ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
        !
        IF (etf_mem) THEN
           etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk, 1)
           etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq, 1)
        ELSE
           nrec = (iq-1) * nksf + ikk
           CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
           nrec = (iq-1) * nksf + ikq
           CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        ENDIF
        !
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
        !
   do iefs = 1, nefsweep            
                !! loop over efs       
        ef1 = efsweep(iefs) !! override for sweeping ef       
        IF ( ( minval ( abs(etf (:, ikk) - ef1) ) .lt. fsthick ) .and. &
             ( minval ( abs(etf (:, ikq) - ef1) ) .lt. fsthick ) ) THEN
           ! instead of using fsthick', fsthick is used here, since only those around the specified Fermi level ef1 will be considered
           !
           fermicount = fermicount + 1
          do ieptemp = 1, neptemp
                eptemp0 = eptemp(ieptemp)
           DO imode = 1, nmodes
              !
              ! the phonon frequency and Bose occupation
              wq = wf (imode, iq)
              wgq = wgauss( -wq/eptemp0, -99)
              wgq = wgq / ( one - two * wgq )
              !
              !  we read the e-p matrix
              !
              IF (etf_mem) THEN
                 epf(:,:) = epf17 ( ik, 1, :, :, imode)
              ELSE
                 nrec = (iq-1) * nmodes * nksqf + (imode-1) * nksqf + ik
                 CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
              ENDIF
              !
              DO ibnd = 1, ibndmax-ibndmin+1
                 !
                 !  the energy of the electron at k (relative to Ef)
                 do ideg = 1, ndeg                       
                         ef0(ideg) = efsweep(iefs) !! override for sweeping ef                   
                         ekk = etf (ibndmin-1+ibnd, ikk) - ef0(ideg)                                                            
                         !
                         DO jbnd = 1, ibndmax-ibndmin+1
                            !
                            !  the fermi occupation for k+q
                            
                                    ekq = etf (ibndmin-1+jbnd, ikq) - ef0(ideg)
                                    wgkq = wgauss( -ekq/eptemp0, -99)  
                                    !
                                    ! here we take into account the zero-point sqrt(hbar/2M\omega)
                                    ! with hbar = 1 and M already contained in the eigenmodes
                                    ! g2 is Ry^2, wkf must already account for the spin factor
                                    !
                                    IF (wq .gt. eps_acustic) THEN
                                        if ( constepm .eqv. .true. ) then
                                                g2 = 0.00000004 / (two*wq) !! assuming g=2.72 meV
                                        else
                                        g2 = abs(epf (jbnd, ibnd))*abs(epf (jbnd, ibnd)) / ( two * wq )
                                      endif
                                    ELSE
                                       g2 = 0.d0
                                    ENDIF
                                    !
                                    ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
                                    ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book 
                                    ! (Many-Particle Physics, 3rd edition)
                                    ! 
                                    weight = wqf(iq) * real (                                        &
                                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw(ideg) )  +  &
                                           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw(ideg) ) ) )
                !                   ecutse needs to be defined if it's used 
                ! @                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                                    !
                                    sigmar(ibndmin-1+ibnd,ikk, ideg, ieptemp, iefs) = sigmar(ibndmin-1+ibnd,ikk, &
                                                                     ideg, ieptemp, iefs) + g2 * weight                                 
                                    !
                                    weight = wqf(iq) * aimag (                                        &
                                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw(ideg) )  +  &
                                           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw(ideg) ) ) ) 
                                    weight_abs = wqf(iq) * aimag (                                        &
                                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw(ideg) )  ) ) 
                                    weight_emi = wqf(iq) * aimag (                                        &
                                         ( ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw(ideg) ) ) ) 
                ! @                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                                    !
                                    sigmai(ibndmin-1+ibnd,ikk, ideg, ieptemp, iefs) = sigmai(ibndmin-1+ibnd,ikk, &
                                                ideg, ieptemp, iefs) + g2 * weight
                                    sigmai_mode_abs(ibndmin-1+ibnd,ikk, ideg, ieptemp, iefs, imode) = &
                                        sigmai_mode_abs(ibndmin-1+ibnd,ikk, ideg, ieptemp, iefs, imode) &
                                                       + g2 * weight_abs
                                    sigmai_mode_emi(ibndmin-1+ibnd,ikk, ideg, ieptemp, iefs, imode) = &
                                        sigmai_mode_emi(ibndmin-1+ibnd,ikk, ideg, ieptemp, iefs, imode) &
                                                        + g2 * weight_emi
                                    !
                                    ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                                    !
                                    weight = wqf(iq) * &
                                         ( (       wgkq + wgq ) * ( (ekk - ( ekq - wq ))**two - degaussw(ideg)**two ) /       &
                                                                  ( (ekk - ( ekq - wq ))**two + degaussw(ideg)**two )**two +  &
                                           ( one - wgkq + wgq ) * ( (ekk - ( ekq + wq ))**two - degaussw(ideg)**two ) /       &
                                                                  ( (ekk - ( ekq + wq ))**two + degaussw(ideg)**two )**two )  
                ! @                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                                    !
                                    zi(ibndmin-1+ibnd,ikk, ideg, ieptemp, iefs) = &
                                         zi(ibndmin-1+ibnd,ikk, ideg, ieptemp, iefs) + g2 * weight
                                 
                            !
                         ENDDO !jbnd
                 enddo !! ndeg
                 !
              ENDDO !ibnd
              !
           ENDDO !imode
          enddo !! ieptemp
           !
        ENDIF ! endif  fsthick
        !
   enddo !! iefs
!!     ENDDO  ! iq's
     !
  ENDDO ! end loop on k
  !
  ! The k points are distributed among pools: here we collect them
  !
  nksqtotf = nkstotf/2 ! odd-even for k,k+q
  !
  ALLOCATE ( etf_all    ( nbndsub, nkstotf ), &
             sigmar_all ( nbndsub, nkstotf, ndeg, neptemp, nefsweep ), &  
             sigmai_all ( nbndsub, nkstotf, ndeg, neptemp, nefsweep ), &  
             zi_all     ( nbndsub, nkstotf, ndeg, neptemp, nefsweep ), &
             sigmai_all_mode_abs ( nbndsub, nkstotf, ndeg, neptemp, nefsweep, nmodes ), &
             sigmai_all_mode_emi ( nbndsub, nkstotf, ndeg, neptemp, nefsweep, nmodes )  )
  IF (.not. ALLOCATED (xkf_all)) then           
    ALLOCATE( xkf_all    ( 3,       nkstotf ))
  ENDIF
  !
  IF (.not. ALLOCATED (sigmai_all_q)) then
    ALLOCATE(sigmai_all_q(nbndsub, nkstotf, ndeg, neptemp, nefsweep)) !! for the sum of selfen_elec contribution from each iq
    sigmai_all_q = zero
  endif
  IF (.not. ALLOCATED (sigmai_all_gpq_abs)) then
    ALLOCATE(sigmai_all_gpq_abs(nbndsub, nkstotf, ndeg, neptemp, nefsweep, nqgroup, nmodes )) !! for the sum of selfen_elec contribution from each iq
    sigmai_all_gpq_abs = zero
  endif
  IF (.not. ALLOCATED (sigmai_all_gpq_emi)) then
    ALLOCATE(sigmai_all_gpq_emi(nbndsub, nkstotf, ndeg, neptemp, nefsweep, nqgroup, nmodes )) !! for the sum of selfen_elec contribution from each iq
    sigmai_all_gpq_emi = zero
  endif
  IF (.not. ALLOCATED (ekk_all)) then
    ALLOCATE(ekk_all(nbndsub, nkstotf)) !! for the electron band energies
    ekk_all = zero
  endif
  
#ifdef __PARA
  !
  ! note that poolgather3 works with the doubled grid (k and k+q)
  ! therefore we need to use the sigma array with both grids, even
  ! though one of them is useless. This should be fixed by modifying
  ! poolgather3 (it's a waste of memory).
  !
  !! potential problem here with multiple deg
  
  CALL poolgather2 ( 3,       nkstotf, nksf, xkf,    xkf_all  )
  CALL poolgather2 ( nbndsub, nkstotf, nksf, etf,    etf_all  )
 ! CALL poolgather3 ( nbndsub, nkstotf, nksf, ndeg, neptemp, nefsweep, sigmar, sigmar_all)
 ! CALL poolgather3 ( nbndsub, nkstotf, nksf, ndeg, neptemp, nefsweep, sigmai, sigmai_all)
 ! CALL poolgather3 ( nbndsub, nkstotf, nksf, ndeg, neptemp, nefsweep, zi,     zi_all)
 ! CALL poolgather4 ( nbndsub, nkstotf, nksf, ndeg, neptemp, nefsweep, nmodes, sigmai_mode_abs, sigmai_all_mode_abs)
 ! CALL poolgather4 ( nbndsub, nkstotf, nksf, ndeg, neptemp, nefsweep, nmodes, sigmai_mode_emi, sigmai_all_mode_emi)
  CALL mp_sum(fermicount, inter_pool_comm)
  !
  ! test output from each pool
  ! DO ik = 1, nksqf
  !    IF (lgamma) THEN
  !       ikk = ik
  !    ELSE
  !       ikk = 2 * ik - 1
  !    ENDIF
  !    WRITE(1000+my_pool_id,'(/5x,"ik = ",i5," coord.: ", 3f9.5)') ik, xkf(:,ikk)
  !    WRITE(1000+my_pool_id,'(6(2x,f12.6))') ( ryd2mev*sigmar(ibnd,ikk), ibnd=ibndmin,ibndmax )
  ! ENDDO
  ! CALL mp_barrier()
  !
#else
  !
  xkf_all = xkf
  etf_all = etf
  sigmar_all = sigmar
  sigmai_all = sigmai
  zi_all     = zi
  sigmai_all_mode_abs = sigmai_mode_abs
  sigmai_all_mode_emi = sigmai_mode_emi
  !
#endif
  !
  !! WRITE(6,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")') 
  !
  !! WRITE( stdout, '(/5x,a,i5,a,i5/)' ) &
  !!  'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nksqtotf*nxqf
  !
    
  DO ik = 1, nksqtotf
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     !
     !! WRITE(stdout,'(/5x,"ik = ",i5," coord.: ", 3f9.5)') ik, xkf_all(:,ikk)
     !! WRITE(stdout,'(5x,a)') repeat('-',67)
     !
     DO ibnd = ibndmin, ibndmax
        !
        ! note that ekk does not depend on q 
        do ideg = 1, ndeg
                do ieptemp =1, neptemp
                        do iefs = 1, nefsweep
                                ekk = etf_all (ibnd, ikk) !! - ef0(ideg)
                                !! with sweeping of multiple efs, it's embiguous to output ekk wrt ef since ef is changing.
                                !
                                ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
                                ! zi_all (ibnd,ikk, ideg) = one / ( one + zi_all (ibnd,ikk, ideg) )
                                !
                        !        WRITE(stdout, 102) ibnd, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ikk), &
                        !            ryd2mev * sigmai_all (ibnd,ikk),  zi_all (ibnd,ikk)
                        !        WRITE(stdout, 103) ik, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ikk), &
                        !              ryd2mev * sigmai_all (ibnd,ikk),  zi_all (ibnd,ikk)
                                ekk_all (ibnd,ikk) = ekk
                                sigmai_all_q (ibnd,ikk, ideg, ieptemp, iefs)  = &
                                        sigmai_all_q (ibnd,ikk, ideg, ieptemp, iefs) &
                                        + sigmai_all (ibnd,ikk, ideg, ieptemp, iefs)                               
                       enddo
                     enddo
                !! sum contribution from each iq to global sigmai
              enddo
       !
      if ( groupq .eqv. .true. ) then
                        !! group phonons into group-ids
                                do ideg = 1, ndeg
                                        do ieptemp =1, neptemp
                                                do iefs = 1, nefsweep
                                                sigmai_all_gpq_abs (ibnd, ikk, ideg, ieptemp, iefs, &
                                                        gpid(iq),:) &
                                                = sigmai_all_gpq_abs (ibnd, ikk, ideg, ieptemp, iefs, gpid(iq), :) &
                                                        + sigmai_all_mode_abs (ibnd,ikk, ideg, ieptemp, iefs, :)
                                                
                                                sigmai_all_gpq_emi (ibnd, ikk, ideg, ieptemp, iefs, &
                                                        gpid(iq),:) &
                                                = sigmai_all_gpq_emi (ibnd, ikk, ideg, ieptemp, iefs, gpid(iq), :) &
                                                        + sigmai_all_mode_emi (ibnd,ikk, ideg, ieptemp, iefs, :)
                                                        !! last column is the nmodes
                                                enddo
                                        enddo
                                enddo
                        else !! no group, default ngroup = 1
                                do ideg = 1, ndeg
                                        do ieptemp =1, neptemp
                                                do iefs = 1, nefsweep
                                                        sigmai_all_gpq_abs (ibnd, ikk, ideg, ieptemp, iefs, 1, :) &
                                                        = sigmai_all_gpq_abs (ibnd, ikk, ideg, ieptemp, iefs, 1, :)&
                                                         + sigmai_all_mode_abs (ibnd,ikk, ideg, ieptemp, iefs, :) !! last column is nmodes
                                                        sigmai_all_gpq_emi (ibnd, ikk, ideg, ieptemp, iefs, 1, :) & 
                                                        = sigmai_all_gpq_emi (ibnd, ikk, ideg, ieptemp, iefs, 1, :)&
                                                         + sigmai_all_mode_emi (ibnd,ikk, ideg, ieptemp, iefs, :) !! single group
                                                        !! last column is the nmodes
                                                enddo
                                        enddo
                                enddo
                        endif
    ENDDO
    
    if ( perphon .eqv. .true. ) then
            do ideg = 1, ndeg
                do ieptemp =1, neptemp
                        do iefs = 1, nefsweep
                                write(72, 107) iq, ik, ideg, ieptemp, iefs, &
                                (ryd2mev *sigmai_all (ibnd,ikk, ideg, ieptemp, iefs) , ibnd=ibndmin,ibndmax)
                        enddo
                enddo
            enddo
          endif
                                
    !! WRITE(stdout,'(5x,a/)') repeat('-',67)
    !
  ENDDO
  !
100 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 FORMAT(5x,'E( ',i3,' )=',f9.3,' eV   Re[Sigma]=',f9.3,' meV   Im[Sigma]=',f9.3,' meV     Z=',f9.3)
103 format(5x,'k( ',i6,' )=',f10.4,' eV   Re[Sigma]=',f10.4,' meV   Im[Sigma]=',f10.4,' meV     Z=',f9.3)
107 format(2i10,3i4, 100e12.4)
  !
  RETURN
  !
  END SUBROUTINE selfen_elec_fly
  
  
  
                                                                           
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_elec !! don't have time to play with this. Will always go with on-the-fly calculations
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the electron self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the electron linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !
  !-----------------------------------------------------------------------
! #include "f_defs.h"
!   USE kinds, ONLY : DP
!   USE cell_base, ONLY : at, bg
!   USE io_global, ONLY : stdout
!   USE phcom, ONLY : lgamma, nmodes
!   USE epwcom, ONLY : nbndsub, lrepmatf, iunepmatf, &
!       fsthick, eptemp, ngaussw, degaussw, iuetf,   &
!       nbndskip, ecutse, parallel_k, ndeg, &
!       parallel_q, epf_mem, etf_mem, eig_read, eps_acustic
!   USE pwcom, ONLY : nelec, ef, isk
!   USE el_phon, ONLY : etf, ibndmin, ibndmax, nksf, etfq, &
!       epf17, wkf, nksqf, nxqf, wf, wqf, xkf, nkstotf, xkf_all
! #ifdef __PARA
!   USE mp,        ONLY : mp_barrier, mp_sum
!   USE mp_global, ONLY : me_pool, inter_pool_comm, my_pool_id
! #endif
!   implicit none
!   !
!   REAL(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
!                               two = 2.d0, zero = 0.d0, pi = 3.14159265358979
!   complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0)
!   integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount, ideg
!   complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
!   REAL(kind=DP) :: g2, ekk, ekq, wq, ef0(ndeg), wgq, wgkq,  &
!        weight, wgauss, dosef(ndeg), dos_ef, sigmar(nbndsub, nksf, ndeg), &
!        sigmai(nbndsub, nksf, ndeg), zi(nbndsub, nksf, ndeg), eptemp0
!   logical :: already_skipped
!   REAL(kind=DP), external :: efermig
!   !
!   ! variables for collecting data from all pools in parallel case 
!   !
!   integer :: nksqtotf
!   REAL(kind=DP), allocatable :: etf_all(:,:), &
!                                 sigmar_all (:,:,:), sigmai_all (:,:,:), zi_all (:,:,:)
!   !
!   WRITE(6,'(/5x,a)') repeat('=',67)
!   WRITE(6,'(5x,"Electron (Imaginary) Self-Energy in the Migdal Approximation")') 
!   WRITE(6,'(5x,a/)') repeat('=',67)
!   !
!   ! loop over temperatures can be introduced
!   !
!   eptemp0 = eptemp(1)
!   !
!   IF ( fsthick .lt. 1.d3 ) &
!     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
!       'Fermi Surface thickness = ', fsthick, ' Ry'
!   WRITE(stdout, '(/5x,a,e18.9,a)' ) &
!     'Golden Rule strictly enforced with T = ',eptemp0, ' Ry'
!   !
!   ! here we take into account that we may skip bands when we wannierize
!   ! (spin-unpolarized)
!   already_skipped = .false.
!   IF ( nbndskip .gt. 0 ) THEN
!      IF ( .not. already_skipped ) THEN 
!         nelec = nelec - two * nbndskip
!         already_skipped = .true.
!         WRITE(stdout, '(/5x,"Skipping the first ",i4," bands:")' ) nbndskip
!         WRITE(stdout, '(/5x,"The Fermi level will be determined with ",f9.5," electrons")' ) nelec
!      ENDIF
!   ENDIF
!   !
!   ! Fermi level and corresponding DOS
!   !
!   do ideg = 1, ndeg
!       ef0(ideg) = efermig(etf,nbndsub,nksf,nelec,wkf,degaussw(ideg),ngaussw,0,isk)  
!         !
!         !   if 'fine' Fermi level differs by more than 250 meV, there is probably
!         !   something wrong with the wannier functions
!         IF ( abs(ef0(ideg) - ef) * ryd2eV .gt. 0.5 .and. (.not.eig_read) ) &
!              CALL errore ('selfen_elec', 'Something wrong with Ef, check MLWFs', -1)
!         !
!         dosef(ideg) = dos_ef (ngaussw, degaussw(ideg), ef0(ideg), etf, wkf, nksf, nbndsub)
!         !   N(Ef) in the equation for lambda is the DOS per spin
!         dosef(ideg) = dosef(ideg) / two
!       
!       !
!         WRITE (6, 100) degaussw(ideg), ngaussw
!         WRITE (6, 101) dosef(ideg), ef0(ideg) * ryd2ev
!         WRITE (6, 101) dosef(ideg), ef  * ryd2ev
!         WRITE (6,'(a)') ' '
!   enddo
!   !
!   sigmar = zero
!   sigmai = zero
!   zi = zero
!   !
!   ! loop over all k points of the fine mesh
!   !
!   fermicount = 0 
!   DO ik = 1, nksqf
!      !
!      IF (lgamma) THEN
!         ikk = ik
!         ikq = ik
!      ELSE
!         ikk = 2 * ik - 1
!         ikq = ikk + 1
!      ENDIF
!      !
!      ! loop over the q-points
!      !
!      DO iq = 1, nxqf
!         !
!         ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
!         !
!         IF (etf_mem) THEN
!            etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk, iq)
!            etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq, iq)
!         ELSE
!            nrec = (iq-1) * nksf + ikk
!            CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
!            nrec = (iq-1) * nksf + ikq
!            CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
!         ENDIF
!         !
!         ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
!         !
!         IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
!              ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
!            !
!            fermicount = fermicount + 1
!            DO imode = 1, nmodes
!               !
!               ! the phonon frequency and Bose occupation
!               wq = wf (imode, iq)
!               wgq = wgauss( -wq/eptemp0, -99)
!               wgq = wgq / ( one - two * wgq )
!               !
!               !  we read the e-p matrix
!               !
!               IF (etf_mem) THEN
!                  epf(:,:) = epf17 ( ik, iq, :, :, imode)
!               ELSE
!                  nrec = (iq-1) * nmodes * nksqf + (imode-1) * nksqf + ik
!                  CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
!               ENDIF
!               !
!               DO ibnd = 1, ibndmax-ibndmin+1
!                  !
!                  !  the energy of the electron at k (relative to Ef)
!                  do ideg = 1, ndeg
!                        ekk = etf (ibndmin-1+ibnd, ikk) - ef0(ideg)
!                        !
!                               DO jbnd = 1, ibndmax-ibndmin+1
!                     !
!                     !  the fermi occupation for k+q
!                     
!                           ekq = etf (ibndmin-1+jbnd, ikq) - ef0(ideg)
!                           wgkq = wgauss( -ekq/eptemp0, -99)  
!                           !
!                           ! here we take into account the zero-point sqrt(hbar/2M\omega)
!                           ! with hbar = 1 and M already contained in the eigenmodes
!                           ! g2 is Ry^2, wkf must already account for the spin factor
!                           !
!                           IF (wq .gt. eps_acustic) THEN
!                              g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
!                           ELSE
!                              g2 = 0.d0
!                           ENDIF
!                           !
!                           ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
!                           ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book 
!                           ! (Many-Particle Physics, 3rd edition)
!                           ! 
!                           weight = wqf(iq) * real (                                        &
!                                ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw(ideg) )  +  &
!                                  ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw(ideg) ) ) )
!       !                   ecutse needs to be defined if it's used 
!       !@                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
!                           !
!                           sigmar(ibndmin-1+ibnd,ikk, ideg) = sigmar(ibndmin-1+ibnd,ikk, ideg) + g2 * weight
!                           !
!                           weight = wqf(iq) * aimag (                                        &
!                                ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw(ideg) )  +  &
!                                  ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw(ideg) ) ) ) 
!       !@                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
!                           !
!                           sigmai(ibndmin-1+ibnd,ikk, ideg) = sigmai(ibndmin-1+ibnd,ikk, ideg) + g2 * weight
!                           !
!                           ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
!                           !
!                           weight = wqf(iq) * &
!                                ( (       wgkq + wgq ) * ( (ekk - ( ekq - wq ))**two - degaussw(ideg)**two ) /       &
!                                                         ( (ekk - ( ekq - wq ))**two + degaussw(ideg)**two )**two +  &
!                                  ( one - wgkq + wgq ) * ( (ekk - ( ekq + wq ))**two - degaussw(ideg)**two ) /       &
!                                                         ( (ekk - ( ekq + wq ))**two + degaussw(ideg)**two )**two )  
!       !@                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
!                           !
!                           zi(ibndmin-1+ibnd,ikk, ideg) = zi(ibndmin-1+ibnd,ikk, ideg) + g2 * weight
! 
!                           !
!                        ENDDO !jbnd
!                      enddo !! ndeg
!                  !
!               ENDDO !ibnd
!               !
!            ENDDO !imode
!            !
!         ENDIF ! endif  fsthick
!         !
!      ENDDO  ! iq's
!      !
!   ENDDO ! end loop on k
!   !
!   ! The k points are distributed among pools: here we collect them
!   !
!   nksqtotf = nkstotf/2 ! odd-even for k,k+q
!   !
!   ALLOCATE ( etf_all    ( nbndsub, nkstotf ), &
!              sigmar_all ( nbndsub, nkstotf, ndeg, neptemp, nefsweep ), &  
!              sigmai_all ( nbndsub, nkstotf, ndeg, neptemp, nefsweep ), &  
!              zi_all     ( nbndsub, nkstotf, ndeg, neptemp, nefsweep )  )
!   IF (.not. ALLOCATED (xkf_all)) then           
!     ALLOCATE( xkf_all    ( 3,       nkstotf ))
!   ENDIF
!              
!   !
! #ifdef __PARA
!   !
!   ! note that poolgather3 works with the doubled grid (k and k+q)
!   ! therefore we need to use the sigma array with both grids, even
!   ! though one of them is useless. This should be fixed by modifying
!   ! poolgather3 (it's a waste of memory).
!   !
!   CALL poolgather2 ( 3,       nkstotf, nksf, xkf,    xkf_all  )
!   CALL poolgather2 ( nbndsub, nkstotf, nksf, etf,    etf_all  )
!   CALL poolgather3 ( nbndsub, nkstotf, nksf, ndeg, neptemp, nefsweep, sigmar, sigmar_all)
!   CALL poolgather3 ( nbndsub, nkstotf, nksf, ndeg, neptemp, nefsweep, sigmai, sigmai_all)
!   CALL poolgather3 ( nbndsub, nkstotf, nksf, ndeg, neptemp, nefsweep, zi,     zi_all)
!   CALL mp_sum(fermicount, inter_pool_comm)
!   !
!   ! test output from each pool
!   ! DO ik = 1, nksqf
!   !    IF (lgamma) THEN
!   !       ikk = ik
!   !    ELSE
!   !       ikk = 2 * ik - 1
!   !    ENDIF
!   !    WRITE(1000+my_pool_id,'(/5x,"ik = ",i5," coord.: ", 3f9.5)') ik, xkf(:,ikk)
!   !    WRITE(1000+my_pool_id,'(6(2x,f12.6))') ( ryd2mev*sigmar(ibnd,ikk), ibnd=ibndmin,ibndmax )
!   ! ENDDO
!   ! CALL mp_barrier()
!   !
! #else
!   !
!   xkf_all = xkf
!   etf_all = etf
!   sigmar_all = sigmar
!   sigmai_all = sigmai
!   zi_all     = zi
!   !
! #endif
!   !
!   WRITE(6,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")') 
!   !
!   WRITE( stdout, '(/5x,a,i5,a,i5/)' ) &
!     'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nksqtotf*nxqf
!   !
!   DO ik = 1, nksqtotf
!      !
!      IF (lgamma) THEN
!         ikk = ik
!         ikq = ik
!      ELSE
!         ikk = 2 * ik - 1
!         ikq = ikk + 1
!      ENDIF
!      !
!      WRITE(stdout,'(/5x,"ik = ",i5," coord.: ", 3f9.5)') ik, xkf_all(:,ikk)
!      WRITE(stdout,'(5x,a)') repeat('-',67)
!      !
!      do ideg = 1, ndeg
!            DO ibnd = ibndmin, ibndmax
!               !
!               ! note that ekk does not depend on q
!               ekk = etf_all (ibnd, ikk) - ef0(ideg)
!               !
!               ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
!               zi_all (ibnd,ikk, ideg) = one / ( one + zi_all (ibnd,ikk, ideg) )
!               !
!       !        WRITE(stdout, 102) ibnd, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ikk), &
!       !            ryd2mev * sigmai_all (ibnd,ikk),  zi_all (ibnd,ikk)
!               WRITE(stdout, 103) ik, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ikk, ideg), &
!                     ryd2mev * sigmai_all (ibnd,ikk, ideg),  zi_all (ibnd,ikk, ideg)
!              !
!           ENDDO
!           WRITE(stdout,'(5x,a/)') repeat('-',67)
!       enddo
!     !
!   ENDDO
!   !
! 100 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
! 101 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
! 102 FORMAT(5x,'E( ',i3,' )=',f9.3,' eV   Re[Sigma]=',f9.3,' meV   Im[Sigma]=',f9.3,' meV     Z=',f9.3)
! 103 format(5x,'k( ',i6,' )=',f10.4,' eV   Re[Sigma]=',f10.4,' meV   Im[Sigma]=',f10.4,' meV     Z=',f9.3)
  !
  RETURN
  !
  END SUBROUTINE selfen_elec


