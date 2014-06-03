 !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine selfen_phon
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the phonon self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the phonon linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  USE phcom,     ONLY : nmodes
  USE epwcom,    ONLY : nbndsub, lrepmatf, iunepmatf, fsthick, &
                        eptemp, ngaussw, degaussw, iuetf,     &
                        wmin, wmax, nw, nbndskip, a2f, epf_mem, etf_mem, &
                        nsmear, delta_smear, eig_read, eps_acustic, egap,&
                        egap0
  USE pwcom,     ONLY : nelec, ef, isk, nbnd
  USE el_phon,   ONLY : epf17, ibndmax, ibndmin, etf, &
                        etfq, wkf, xqf, wqf, nksf, nxqf,   &
                        nksqf, wf, nkstotf, xkf, xqf, dmef, &
                        lambda_all, lambda_v_all
#ifdef __PARA
  USE mp,        ONLY : mp_barrier,mp_sum
  USE mp_global, ONLY : me_pool,inter_pool_comm,my_pool_id
#endif
  !
  implicit none
  !
  real(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0), czero=(0.d0,0.d0)
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount, ismear
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgkk, wgkq, lambda, lambda_v, &
                   weight, wgauss, dosef, dos_ef, w0g1, w0g2, w0gauss,    &
                   gamma(nmodes),gamma_v(nmodes), degaussw0, eptemp0, lambda_tot, degausswi
  !
  real(kind=DP), external :: efermig
  real(kind=DP) :: gamma_all  ( nmodes, nxqf, 10 )
  real(kind=DP) :: gamma_all_v  ( nmodes, nxqf, 10 )
  real(kind=DP) :: coskkq(nbndsub, nbndsub)
  real(kind=DP) :: DDOT,  vkk(3,nbndsub), vkq(3,nbndsub)
  !
  !
  !
  degausswi = degaussw(1) !! Do not loop over degauss for phonon at this time
  
  WRITE(6,'(/5x,a)') repeat('=',67)
  WRITE(6,'(5x,"Phonon (Imaginary) Self-Energy in the Migdal Approximation")') 
  WRITE(6,'(5x,a/)') repeat('=',67)
  !
  IF ( fsthick .lt. 1.d3 ) &
     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
     'Fermi Surface thickness = ', fsthick, ' Ry'
  WRITE(stdout, '(/5x,a,e18.9,a)' ) &
     'Golden Rule strictly enforced with T = ',eptemp(1), ' Ry'
  ! here we take into account that we may skip bands when we wannierize
  ! (spin-unpolarized)
  ! 
  IF (nbndskip.gt.0) then
     nelec = nelec - two * nbndskip
     WRITE(stdout, '(/5x,"Skipping the first ",i4," bands:")' ) nbndskip
     WRITE(stdout, '(/5x,"The Fermi level will be determined with ",f9.5," electrons")' ) nelec     
  ENDIF
  !
  ! here we loop on smearing values JN - this may be an option later on 
  !
  DO ismear = 1, nsmear
     !
     degaussw0 = (ismear-1) * delta_smear + degausswi
     eptemp0 = (ismear-1) * delta_smear + eptemp(1)
     !
     IF (.not.ALLOCATED (lambda_all) .and. .not.ALLOCATED(lambda_v_all)) then
        ALLOCATE(lambda_all(nmodes, nxqf, nsmear))
        ALLOCATE(lambda_v_all(nmodes, nxqf, nsmear))
     ENDIF
     gamma_all =  0.d0
     gamma_all_v = 0.d0
     lambda_all = 0.d0
     lambda_v_all = 0.d0
     ! 
     !
     ! Fermi level and corresponding DOS
     !
     !   Note that the weights of k+q points must be set to zero here
     !   no spin-polarized calculation here
     ef0 = efermig(etf,nbndsub,nksf,nelec,wkf,degaussw0,ngaussw,0,isk)
     !
     !   if 'fine' Fermi level differs by more than 250 meV, there is probably
     !   something wrong with the wannier functions
     IF (abs(ef0 - ef) * ryd2eV .gt. 0.5 .and. (.not.eig_read) ) &
        CALL errore ('selfen_phon', 'Something wrong with Ef, check MLWFs', 1)
     !
     dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub)
     !   N(Ef) in the equation for lambda is the DOS per spin
     dosef = dosef / two
     !
     WRITE (6, 100) degaussw0, ngaussw
     WRITE (6, 101) dosef, ef0 * ryd2ev
     !
     ! loop over all q points of the fine mesh (this is in the k-para case 
     ! it should always be k-para for selfen_phon)
     ! 
     DO iq = 1, nxqf
        !
        CALL start_clock('PH SELF-ENERGY')
        !
        fermicount = 0
        gamma = zero
        gamma_v = zero
        !
        DO ik = 1, nksqf !nksqf: number of kpoint blocks in the pool
           !
           ikk = 2 * ik - 1
           ikq = ikk + 1
           !
           coskkq = 0.d0
           DO ibnd = 1, nbndsub
              DO jbnd = 1, nbndsub
                 ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
                 ! 1/m  = 2 in Rydberg atomic units
                 vkk(:, ibnd ) = 2.0 * REAL (dmef (:, ibnd, ibnd, ikk ) )
                 vkq(:, jbnd ) = 2.0 * REAL (dmef (:, jbnd, jbnd, ikq ) )
                 IF ( abs ( DDOT(3,vkk(:,ibnd), 1, vkk(:,ibnd), 1) ) .gt. 1.d-4 ) &
                    coskkq(ibnd, jbnd ) = DDOT(3, vkk(:,ibnd ), 1, vkq(:,jbnd),1) / &
                                          DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd),1)
              ENDDO
           ENDDO
           !
           ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
           !
           IF (etf_mem) then
              etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk, iq)
              etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq, iq)
           ELSE
              nrec = (iq-1) * nksf + ikk
              CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
              nrec = (iq-1) * nksf + ikq
              CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
           ENDIF
           !
           ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
           IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
                ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) then
              !
              fermicount = fermicount + 1
              !
              DO imode = 1, nmodes
                 !
                 ! the phonon frequency
                 wq = wf (imode, iq)
                 !
                 !  we read the e-p matrix from disk / memory
                 !
                 IF (etf_mem) then
                    epf(:,:) = epf17 ( ik, iq, :, :, imode)
                 ELSE
                    nrec = (iq-1) * nmodes * nksqf + (imode-1) * nksqf + ik
                    CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
                 ENDIF
                 !
                 DO ibnd = 1, ibndmax-ibndmin+1
                    !
                    !  the fermi occupation for k
                    ekk = etf (ibndmin-1+ibnd, ikk) - ef0
                    wgkk = wgauss( -ekk/eptemp0, -99)
                    !
                    DO jbnd = 1, ibndmax-ibndmin+1
                       !
                       !  the fermi occupation for k+q
                       ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                       wgkq = wgauss( -ekq/eptemp0, -99)  
                       !
                       ! here we take into account the zero-point sqrt(hbar/2M\omega)
                       ! with hbar = 1 and M already contained in the eigenmodes
                       ! g2 is Ry^2, wkf must already account for the spin factor
                       !
                       ! NON FUNZIONA SCAMBIANDO i,j
                       ! the coupling from Gamma acoustic phonons is negligible
                       IF ( wq .gt. eps_acustic ) THEN
                          g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
                       ELSE
                          g2 = 0.d0
                       ENDIF
                       !
                       ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
                       ! This is the imaginary part of the phonon self-energy, sans the matrix elements
                       
                       weight = wkf (ikk) * (wgkk - wgkq) * &
                          aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) ) 
                       
                       ! the below expression is positive-definite, but also an approximation
                       ! which neglects some fine features
                       !
                      ! w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
                      ! w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
                      ! weight = pi * wq * wkf (ikk) * w0g1 * w0g2
                       !
                       gamma(imode) =   gamma   (imode) + weight * g2 
                       gamma_v(imode) = gamma_v (imode) + weight * g2 * (1-coskkq(ibnd, jbnd) ) 
                       !
                    ENDDO ! jbnd
                 ENDDO   ! ibnd
                 !
              ENDDO ! loop on q-modes
              !
              !
           ENDIF ! endif fsthick
           !
           CALL stop_clock('PH SELF-ENERGY')
           !
        ENDDO ! loop on k
#ifdef __PARA
        !
        ! collect contributions from all pools (sum over k-points)
        ! this finishes the integral over the BZ  (k)
        !
        CALL mp_sum(gamma,inter_pool_comm) 
        CALL mp_sum(gamma_v,inter_pool_comm) 
        CALL mp_sum(fermicount, inter_pool_comm)
        !
#endif
        !
        WRITE(6,'(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') iq, xqf(:,iq) , wqf(iq)
        WRITE(6,'(5x,a)') repeat('-',67)
        lambda_tot = zero
        DO imode = 1, nmodes
           ! 
           wq = wf (imode, iq)
           lambda = zero
           lambda_v = zero
           IF ( sqrt(abs(wq)) .gt. eps_acustic ) lambda   = gamma(imode) / pi / wq**two / dosef
           IF ( sqrt(abs(wq)) .gt. eps_acustic ) lambda_v = gamma_v(imode) / pi / wq**two / dosef
           lambda_tot = lambda_tot + lambda
           !
           gamma_all( imode, iq, ismear ) = gamma(imode)
           gamma_all_v( imode, iq, ismear ) = gamma_v(imode) / pi / wq**two / dosef
           lambda_all( imode, iq, ismear ) = lambda
           lambda_v_all( imode, iq, ismear ) = lambda_v
           !
           WRITE(6, 102) imode, lambda, ryd2mev * gamma(imode), ryd2mev * wq
  call flush(6)
        ENDDO
        !
        WRITE(6,103) lambda_tot
        WRITE(6,'(5x,a/)') repeat('-',67)
        !
        ! test ONLY
#ifdef __PARA
!        if (me.eq.1) & 
        IF (me_pool == 0) &
#endif
        !
        !     
        WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
             'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkstotf/2
        !
     ENDDO ! loop on q
     !
  ENDDO !smears
  ! generate the Eliashberg spectral function
  !
  IF (a2f) call eliashberg_a2f( lambda_all(:,:,1), lambda_v_all(:,:,1))
  !
100 format(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 format(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 format(5x,'lambda( ',i3,' )=',f9.3,'   gamma=',f9.3,' meV','   omega=',f9.3,' meV')
103 format(5x,'lambda( tot )=',f9.3)
  !
  end subroutine selfen_phon
  !
  !-----------------------------------------------------------------------
  subroutine selfen_phon_fly (iq )
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the phonon self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the phonon linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation.  This routine is similar to the one above
  !  but it is ONLY called from within ephwann_shuffle and calculates 
  !  the selfenergy for one phonon at a time.  Much smaller footprint on
  !  the disk
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  use phcom,     ONLY : nmodes
  use epwcom,    ONLY : nbndsub, lrepmatf, iunepmatf, fsthick, &
                        eptemp, ngaussw, degaussw, iuetf,     &
                        wmin, wmax, nw, nbndskip, epf_mem, etf_mem, &
                        ndeg, delta_smear, eig_read, eps_acustic, efsweep, nefsweep,&
                        neptemp, ldoubledelta, egap, egap0, ltetra_phon, &
                        etf_global, lower_bnd, upper_bnd
  use pwcom,     ONLY : nelec, isk, nbnd, ef
  use el_phon,   ONLY : epf17, ibndmax, ibndmin, etf, &
                        etfq, wkf, xqf, wqf, nksf, nxqf,   &
                        nksqf, wf, nkstotf, xkf, xqf, &
                        lambda_all, lambda_v_all, nrr_k, &
                        dmef, ndegen_k, irvec
  USE tetrahedron
#ifdef __PARA
  use mp,        ONLY : mp_barrier,mp_sum, mp_bcast
  use mp_global, ONLY : mpime,me_pool,inter_pool_comm,my_pool_id
  USE io_global, ONLY : ionode_id
#endif
  !
  implicit none
  !
  real(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0), czero=(0.d0,0.d0)
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, ismear, &
              iefs, ieptemp, ierr
  complex(kind=DP) epf
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgkk, wgkq, lambda, lambda_v, &
     weight, wgauss, dosef, dos_ef, w0g1, w0g2, w0gauss,   &
     gamma(nmodes,ndeg,nefsweep,neptemp),gamma_v(nmodes,ndeg,nefsweep,neptemp),&
     degaussw0, eptemp0, lambda_tot
  !
  real(kind=DP), external :: efermig
  real(kind=DP) :: coskkq(nbndsub, nbndsub)
  real(kind=DP) :: DDOT,  vkk(3,nbndsub), vkq(3,nbndsub)
  real(kind=DP) :: deltaegap !The adjusted band gap 
  !
  !local variables for tetrahedron
  real(kind=DP) :: eprime(nkstotf/2)
  real :: fk
  integer :: iktetra
  integer :: iii,jjj !test purpose
  integer :: itetra_(2,30,nksqf) !itetra in this pool

  IF (iq.eq.1) then 
     WRITE(6,'(/5x,a)') repeat('=',67)
     WRITE(6,'(5x,"Phonon (Imaginary) Self-Energy in the Migdal Approximation (on the fly)")') 
     WRITE(6,'(5x,a/)') repeat('=',67)
     !
     IF ( fsthick.lt.1.d3 ) &
          WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Fermi Surface thickness = ', fsthick, ' Ry'
     !WRITE(stdout, '(/5x,a,e18.9,a)' ) &
     !     'Golden Rule strictly enforced with T = ',eptemp(1), ' Ry'
     !
     ! here we take into account that we may skip bands when we wannierize
     ! (spin-unpolarized)
     ! 
     ! could be done for iq = 1
     !IF (.not.ALLOCATED (lambda_all) .and. .not.ALLOCATED(lambda_v_all)) then
     !   ALLOCATE(lambda_all(nmodes, nxqf, ndeg))
     !   ALLOCATE(lambda_v_all(nmodes, nxqf, ndeg))
     !   lambda_all   = zero
     !   lambda_v_all = zero
     !ENDIF
     !
     IF (nbndskip.gt.0) then
        nelec = nelec - two * nbndskip
        WRITE(stdout, '(/5x,"Skipping the first ",i4," bands:")' ) nbndskip
        WRITE(stdout, '(/5x,"The Fermi level will be determined with ",f9.5," electrons")' ) nelec     
     ENDIF
     !
  ENDIF
  !
  ! Adjust the bandgap to the imposed value
  deltaegap=0
  IF (egap0 .ne. 0) THEN
     WRITE(stdout,'(/5x,a)') 'Imposed band gap is used'
     WRITE(stdout,'(/5x,a,f9.3)') 'The calculated band gap is (eV)',egap*ryd2ev
     WRITE(stdout,'(/5x,a,f9.3)') 'The imposed band gap is (eV)',egap0*ryd2ev
     deltaegap=egap0-egap
  ENDIF
  !!!!
!  WRITE(stdout,'(/5x,a)') 'Starting the velocity loop'
  DO ik = 1, nksqf !start the velocity loop
     !
     ikk = 2 * ik - 1
     ikq = ikk + 1
     ! 
     coskkq = 0.d0
     DO ibnd = 1, nbndsub
        DO jbnd = 1, nbndsub
           ! coskkq = (vk dot vkq) / |vk|^2  appears in Grimvall 8.20
           ! this is different from :   coskkq = (vk dot vkq) / |vk||vkq|
           ! In principle the only coskkq contributing to lambda_tr are both near the
           ! Fermi surface and the magnitudes will not differ greatly between vk and vkq
           ! we may implement the approximation to the angle between k and k+q vectors also 
           ! listed in Grimvall
           !
           ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
           vkk(:, ibnd ) = 2.0 * REAL (dmef (:, ibnd, ibnd, ikk ) )
           vkq(:, jbnd ) = 2.0 * REAL (dmef (:, jbnd, jbnd, ikq ) )
           if ( abs ( vkk(1,ibnd)**2 + vkk(2,ibnd)**2 + vkk(3,ibnd)**2) .gt. 1.d-4) &
                coskkq(ibnd, jbnd ) = DDOT(3, vkk(:,ibnd ), 1, vkq(:,jbnd),1)  / &
                DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd),1)
        ENDDO
     ENDDO
     ! when we see references to iq for file readinq, it is always = 1 for on the fly calculations
     IF (etf_mem) then
        etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk,  1)
        etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq,  1)
     ELSE
        nrec = (iq-1) * nksf + ikk
        nrec = ikk
        CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        nrec = (iq-1) * nksf + ikq
        nrec = ikq
        CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
     ENDIF
     !
  ENDDO !End of the velocity loop
  !Velocity loop finished
!  WRITE(stdout,'(/5x,a)') 'Finished the velocity loop'

  CALL start_clock('PH SELF-ENERGY')

  gamma = zero
  gamma_v = zero
  
  IF (ltetra_phon) THEN
      IF (upper_bnd/2-(lower_bnd+1)/2+1 .ne. nksqf) CALL errore('tetra','k-points inconsistent')
      itetra_(:,:,:) = itetra(:,:,(lower_bnd+1)/2:upper_bnd/2)
  ENDIF
!  WRITE(stdout,'(/5x,a)') 'Here reached, 5'
!  WRITE(stdout,'(/5x,a,i8)') 'Total number of kpoints: ', nkstotf

  DO imode = 1, nmodes !loop over modes
     !
     ! the phonon frequency
     wq = wf (imode, iq)
     !
     !
     DO ibnd = 1, ibndmax-ibndmin+1
        DO jbnd = 1, ibndmax-ibndmin+1
           
!           WRITE(stdout,'(/5x,a,i4,a,i4,a,i4,a,i4)') 'Calculating iq =',iq,'imode=', imode, 'm =',ibnd,'n =',jbnd    
           !
           !Here the tetrahedron are updates with eigenevalues of electrons
           !
           IF (ltetra_phon) THEN
                   !CALL mp_barrier()
                   DO ik = 1, nkstotf/2
     !                  eprime(ik)=etf_global(ibndmin+ibnd-1,2*ik-1)-etf_global(ibndmin+jbnd-1,2*ik)
                       eprime(ik)=etf_global(ibndmin+jbnd-1,2*ik)-etf_global(ibndmin+ibnd-1,2*ik-1)
     !                  WRITE(stdout,'(/5x,a,i8)') 'Here reached',ik
                       !IF (ik .le. 100) WRITE(stdout,'(2f13.7)') etf_global(ibndmin+ibnd-1,2*ik-1), etf_global(ibndmin+jbnd-1,2*ik)
                   ENDDO
                 !  CALL mp_barrier()
                   !fill eigenvalues to the tetrahedron
                   CALL eigen_tet(ntetra,eprime,tetra_i,tetra_w,nkstotf/2) 
                 !      WRITE(stdout,'(/5x,a)') 'Here reached,8'
                   !fill weights to the tetrahedron
                   CALL weight_tet(nkstotf/2,ntetra,wq,tetra_i,tetra_w,tetra_c,wkt)
                  !     WRITE(stdout,'(/5x,a)') 'Here reached,9'
!!!!#ifdef __PARA
                 IF (mpime .eq. ionode_id) THEN
                   IF ((iq .eq. 10) .and. (ibnd .eq. 1) .and. (jbnd .eq. 1)) THEN
                       IF (imode .eq. 3) THEN
                           OPEN(UNIT=11111,FILE='tetratest.dat',STATUS='replace')
                           WRITE(11111,'(a,f13.7)') 'Phonon energy = ',wq
                           DO iii=1,ntetra
                               DO jjj=1,4
                                   WRITE(11111,'(i4,i6,i4,i6,2f13.7)') imode,iii,jjj,tetra_i(jjj,iii),tetra_w(jjj,iii),tetra_c(jjj,iii)
                               ENDDO
                           ENDDO
                       ENDIF
                    ENDIF
                ENDIF
!#endif
         ENDIF 
           !
                  ! CALL mp_barrier()
                   !    WRITE(stdout,'(/5x,a)') 'Here reached,10'
           DO ik = 1, nksqf !k loop
              ikk = 2 * ik - 1
              ikq = ikk + 1
              IF ( (abs(etf (ibnd, ikk) - ef) .lt. fsthick ) .and. &
              (abs(etf (jbnd, ikq) - ef) .lt. fsthick ) ) then
              DO ismear = 1, ndeg  !degauss loop
                 degaussw0 = degaussw(ismear)
                 DO iefs = 1, nefsweep !Ef loop
                     ef0 = efsweep(iefs) !!overide for sweeping ef
                     DO ieptemp = 1, neptemp !temperature loop
                        eptemp0 = eptemp(ieptemp)
                        !
                        !  the fermi occupation for k
                        ekk = etf (ibndmin-1+ibnd, ikk) - ef0
                        ! impose the band gap
                        IF (etf(ibndmin-1+ibnd, ikk) .ge. ef) ekk = ekk + deltaegap 
                        wgkk = wgauss( -ekk/eptemp0, -99)
                        !
                        !  the fermi occupation for k+q
                        ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                        ! impose the band gap
                        IF (etf(ibndmin-1+jbnd, ikq) .ge. ef) ekq = ekq + deltaegap 
                        wgkq = wgauss( -ekq/eptemp0, -99)  
                        !
                        ! here we take into account the zero-point sqrt(hbar/2M\omega)
                        ! with hbar = 1 and M already contained in the eigenmodes
                        ! g2 is Ry^2, wkf must already account for the spin factor
                        !
                        ! NON FUNZIONA SCAMBIANDO i,j
                        ! the coupling from Gamma acoustic phonons is negligible
                        IF (etf_mem) then
                           epf = epf17 ( ik,  1, jbnd, ibnd, imode)
                        ELSE !if etf_mem is false, here modification is needed.-BL
                          ! nrec = (imode-1) * nksqf + ik
                          ! CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
                        ENDIF

                        IF ( wq .gt. eps_acustic ) THEN
                           g2 = abs(epf)**two / ( two * wq )
                        ELSE
                           g2 = 0.d0
                        ENDIF
                        !
                        ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
                        ! This is the imaginary part of the phonon self-energy, sans the matrix elements
                        !
                        IF (.not. ldoubledelta .and. .not. ltetra_phon) THEN
                            weight = wkf (ikk) * (wgkk - wgkq) * &
                            aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) ) 
                            gamma(imode,ismear,iefs,ieptemp) = gamma(imode,ismear,iefs,ieptemp) + weight * g2 
                            gamma_v(imode,ismear,iefs,ieptemp) = gamma_v (imode,ismear,iefs,ieptemp) + weight * g2 * (1-coskkq(ibnd, jbnd) ) 
                        ! WRITE(6,'(/5x,a)') 'Exact FD occupation number is used'
                        ELSEIF (.not. ltetra_phon) THEN !double-delta approximation
                        ! the below expression is positive-definite, but also an approximation
                        ! which neglects some fine features
                        !
                            w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
                            w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
                            weight = pi * wq * wkf (ikk) * w0g1 * w0g2
                            gamma(imode,ismear,iefs,ieptemp) = gamma(imode,ismear,iefs,ieptemp) + weight * g2 
                            gamma_v(imode,ismear,iefs,ieptemp) =gamma_v(imode,ismear,iefs,ieptemp) + weight * g2 * (1-coskkq(ibnd, jbnd) ) 
                        !WRITE(6,'(/5x,a)') 'Double delta approximation is used'
                        ELSE !tetrahedron is used
                   !         CALL mp_barrier()
                    !        WRITE(stdout,'(/5x,a)') 'Here reached,11'
                            fk =  g2 * (wgkk - wgkq) !Here the weight of k-point is already taken care of in the weights of tetrahedron. Don't overcound
                            iktetra = 1 !search for all tetrahedron that this k-point belongs to
                            DO WHILE (itetra_(1,iktetra,ik) .ne. 0)
                                gamma(imode,ismear,iefs,ieptemp) =gamma(imode, ismear, iefs, ieptemp) + tetra_c(itetra_(2,iktetra,ik),itetra_(1,iktetra,ik)) * fk * pi 
                                gamma_v(imode,ismear,iefs,ieptemp) =gamma_v (imode, ismear, iefs, ieptemp) + tetra_c(itetra_(2,iktetra,ik),itetra_(1,iktetra,ik)) * fk * pi * (1-coskkq(ibnd, jbnd) ) 
                                iktetra = iktetra + 1 !look for the next tetrahedra
                                if (iktetra .gt. 30) CALL errore('tetra','too many tetrahedron associated with a k-point!',1)
                            ENDDO    
                        ENDIF
                     ENDDO !loop on temperature  
                 ENDDO !loop on Ef  
              ENDDO !loop on degauss     
              ENDIF
           ENDDO !loop on k   
        ENDDO ! jbnd
     ENDDO   ! ibnd
     !
!  WRITE(stdout,'(/5x,a)') 'Here reached, 6'

  ENDDO ! loop on q-modes

CALL stop_clock('PH SELF-ENERGY')

#ifdef __PARA
     !
     ! collect contributions from all pools (sum over k-points)
     ! this finishes the integral over the BZ  (k)
     CALL mp_barrier()
     CALL mp_sum(gamma,inter_pool_comm) 
     CALL mp_sum(gamma_v,inter_pool_comm) 
     !
#endif
     !
  DO ismear=1,ndeg
      DO iefs=1,nefsweep
          dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub) 
          dosef = dosef / two
          DO ieptemp=1,neptemp
              WRITE(6,'(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5, " ieptemp: ", i3, " iefs: ",i3, " ismear: ",i3)') iq, xqf(:,iq) , wqf(iq), ieptemp, iefs, ismear
              WRITE(6,'(5x,a)') repeat('-',67)
!              lambda_tot = 0
              DO imode = 1, nmodes
!              ! 
                  wq = wf (imode, iq)
!                  lambda   = zero
!                  lambda_v = zero
!                  IF ( sqrt(abs(wq)) .gt. eps_acustic ) lambda   = gamma(imode,ismear,iefs,ieptemp)   / pi / wq**two / dosef
!                  IF ( sqrt(abs(wq)) .gt. eps_acustic ) lambda_v = gamma_v(imode) / pi / wq**two / dosef
!                  lambda_tot = lambda_tot + lambda
        !
!                  lambda_all( imode, iq, ismear ) = lambda
!                  lambda_v_all( imode, iq, ismear ) = lambda_v
        !
        !
                  WRITE(6, 102) imode, ryd2mev * gamma(imode,ismear, iefs, ieptemp), ryd2mev * wq
               ENDDO
!               WRITE(6, 103) lambda_tot
               WRITE(6,'(5x,a/)') repeat('-',67)
     !
     !
     ENDDO !ieptemp
   ENDDO ! iefs
  ENDDO !smears
  !
100 format(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 format(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 format(5x,'Mode index=',i3,' gamma=',f13.7,' meV',' omega=',f13.7,' meV')
103 format(5x,'lambda( sum )=',f13.7)
  !
!IF (ltetra_phon) THEN
!   DEALLOCATE(tetra_i,tetra_w, tetra_c,itetra,wkt)
!ENDIF
  return
end subroutine selfen_phon_fly
!

