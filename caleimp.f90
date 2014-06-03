!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE caleimp()
  !-----------------------------------------------------------------------
  !
  ! Electron-impurity calculation from data saved in diffv.dat
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : celldm, omega, ibrav, alat
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, pmass
  USE gvecs, ONLY: doublegrid
  USE fft_base, ONLY : dfftp, dffts
  USE noncollin_module, ONLY : nspin_mag
  USE dynmat, ONLY : dyn, w2
  USE qpoint, ONLY : xq
  USE modes,  ONLY : npert, nirr
  USE control_ph, ONLY : trans
  USE units_ph, ONLY : iudyn, lrdrho, iudvscf
  !
  IMPLICIT NONE
  !
  INTEGER :: irr, imode0, ipert, is, ir, stat, nx, ny, nz, nnrdv, analydV, &
             ix, iy, iz
  real(DP) :: LdV,tmpr, bohr=0.52917720859d0, lambda, pref, rabs(3), rlen
  character :: line*130
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  REAL(DP), ALLOCATABLE :: diffV(:)

	write(*,*) 'Ready to begin compute_Rate()'
 CALL compute_rate ()
	write(*,*) 'Do it in three steps:'
	write(*,*) 'Now the coarse grid elimpmat has been stored in disk'
	write(*,*) 'run matlab in the current folder for 3d spline interpolation'
	write(*,*) 'run another postproc code to use the fit grid elimpmat for rate'
	call exit(1)
     !
!   DEALLOCATE (diffV)

  !
  ! now read the eigenvalues and eigenvectors of the dynamical matrix
  ! calculated in a previous run
!   !
!   IF (.NOT.trans) CALL readmat (iudyn, ibrav, celldm, nat, ntyp, &
!        ityp, omega, pmass, tau, xq, w2, dyn)
  !
  CALL stop_clock ('eimp')
  RETURN
END SUBROUTINE caleimp
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine compute_rate()
USE qpoint,   ONLY : nksq, npwq, ikks
USE kinds,    ONLY : DP
USE klist,    ONLY: xk, wk
USE cell_base, ONLY : alat, at, bg
USE fft_base, ONLY : dffts
USE save_ph,     ONLY : tmp_dir_save
USE io_files,    ONLY : prefix, tmp_dir, iunigk
USE wvfct, ONLY: nbnd, npw, npwx, igk
USE lsda_mod, ONLY: lsda, current_spin, isk
USE noncollin_module, ONLY : noncolin, npol, nspin_lsda, nspin_mag
USE wavefunctions_module,  ONLY: evc
USE units_ph, ONLY : iubar, lrbar, lrwfc, iuwfc
USE control_ph, ONLY : trans, lgamma, tmp_dir_phq, xmldyn
USE mp_global, ONLY: intra_pool_comm
USE mp,        ONLY: mp_sum
USE parameters,  ONLY : npk
USE symm_base,   ONLY: s, irt, nsym, invs
USE io_global,   ONLY : stdout, ionode, ionode_id
implicit none
! LOCAL variables

INTEGER :: iuna2Fsave  = 40

INTEGER ::  nscf
REAL(DP) :: kimkj(3), LdV
real(dp) :: ks, kBT, e0, pi, mm, ee, hbar, kimkjm, ki(3), kj(3), kim, kjm
REAL(DP), PARAMETER :: ryd2mev = 13605.8, bohr=0.52917720859d0
INTEGER :: nrec, ik, ikk, ikq, ipert, mode, ibnd, jbnd, ir, ig, &
ios, jk, jkk, analyelmat
COMPLEX(DP) , ALLOCATABLE :: evir (:,:),evjr (:,:), evrtot (:,:,:,:), evr(:,:), psir_ex(:,:), gf(:,:,:,:), gf1(:, :,:,:, :)
complex(dp), allocatable :: eimpmt(:)
real(dp), allocatable :: eig(:,:)
integer :: parabo, useakkp, bndsta, bndend, i, isig, nsig, flag, bndcbm, kcbm, numnd, ind
real(dp) :: degauss1, deg(8), thrswindow, etk, etkj, w0g1, costheta, wqa, cbm, thrs, pref, wkinp, xksft(3), vthrs
complex(dp) :: elimpplwv, akkp
complex(DP), parameter :: cim = (0.d0,1.d0)
character*130 :: filename, lines
real(dp), allocatable :: nd(:), Ld(:)

INTEGER, ALLOCATABLE :: eqBZ(:), sBZ(:), map(:)
REAL(DP), allocatable :: xkfit(:,:), xkfitt(:,:), vik(:,:,:)
REAL(DP), allocatable, target :: etfit(:,:), wkfit(:)
INTEGER :: nksfit, nk1fit, nk2fit, nk3fit, nkfit
INTEGER, allocatable :: eqkfit(:),sfit(:)
INTEGER :: nks_real, ispin, nksfit_real
integer ::  ikx, iky, ikz, jkx, jky, jkz, nxyz, &
            ikxp1, ikyp1, ikzp1, indxp, indyp, indzp, tmpi, tmpi1, &
            ikxm1, ikym1, ikzm1, indxm, indym, indzm, l
real(dp) :: etknbp(3), etknbm(3), v1(3), v2(3), absdk(3), vec(3,3), vtmp(3)
LOGICAL  :: exst, xmldyn_save
  
  
bndsta = 5
bndend = 6

nscf = floor(dffts%nnr**(1.0/3))
ALLOCATE (evr(dffts%nnr, npol),evir(dffts%nnr, npol),evjr(dffts%nnr, npol), evrtot(nksq, bndend-bndsta + 1,dffts%nnr, npol)) !! npol =1 !!, evrtot(nksq, nbnd,nnrdv, npol)

ks=11.68
kBT = 0.027 !! 300K ~ 0.027 eV
pref = (1d-34)*(1d-34)*(1d10)*(1d10)/2.0/(1.0*9.1d-31)/(1.6d-19)/13.6 !! E= hbar^2 * k^2/ 2m, in Ry, k in Ang^-1
!! in the excel file, it's seen taking m*=1.0 gives best fit to band shape around cbm

vthrs =0.002 !! number obtained by looking at the actual "v" output and I know cbm v should be 0

nsig = 5

deg(1) = 0.0001 !! Ry
deg(2) = 0.0005 !! Ry
deg(3) = 0.001 !! Ry
deg(4) = 0.002 !! Ry
deg(5) = 0.005 !! Ry


if( parabo .eq. 1) bndend = 5

open(79, file="calinput.txt", status="unknown")
read(79, *) parabo !! parabo = 1
read(79, *) useakkp  !! whether to use akkp or just take plane wave solution
read(79, *) bndcbm !! read in the band where cbm is
read(79, *) thrs !! only energy between cbm and cbm+thrs will be computed , in eV
read(79, *) wkinp !! input for weights per k point. Otherwise the one from scf is wrong
read(79, *) numnd !! number of impurity density to be studied

allocate(nd(numnd), Ld(numnd), eimpmt(numnd))
read(79, *) (nd(ind), ind = 1, numnd)  !! impurity density in m^-3
do ind = 1, numnd
	Ld(ind)= sqrt(ks*8.854e-12*kBT/(1.6e-19*nd(ind)))/(1.0e-10) !! take Ld from 300K, based on nd. Debye length Ld in Ang
enddo

write(*,*) 'Read calinput.txt:'
if(parabo .eq. 1) then
	write(*,*) 'Use parabollic bandstr, origin to be shifted to k_cbm'
else
	write(*,*) 'DFT bandstructure is used'
endif
if(useakkp .eq. 1) then
	write(*,*) 'Add Akk_prime coefficient to multiply with planewave solution '
else
	write(*,*) 'Planewave solution only'
endif
write(*,*) 'The', bndcbm,'th band contains cbm (from input)'
write(*,*) 'Only energy within', thrs,' eV from cbm will be output for final rate'
do ind = 1, numnd
	write(*,*) 'Debye length:', Ld(ind),' Ang, impurity density:', nd(ind), 'm^-3'
enddo
write(*,*) 'Weights per kpt: ', wkinp
open(80, file="seltofullk.txt", status = "unknown")
read(80, *)  !! should read in nksq 
read(80, *)
allocate(map(nksq))
do ik = 1, nksq
	read(80, *) tmpi, tmpi1
	map(tmpi) = tmpi1
enddo

! read eigenvalues for the dense grid
  ! FIXME: this might be done from the xml file, not from a specialized file
  ! parallel case: only first node reads
  !
  IF ( ionode ) THEN
     tmp_dir=tmp_dir_save
     write(*,*) 'tmp_dir: ', tmp_dir
     CALL seqopn( iuna2Fsave, 'a2Fsave', 'FORMATTED', exst )
     tmp_dir=tmp_dir_phq
     READ(iuna2Fsave,*) ibnd, nksfit
  END IF
    
  allocate (etfit(nbnd,nksfit), xkfit(3,nksfit), wkfit(nksfit))
  !
  IF ( ionode ) THEN
     READ(iuna2Fsave,*) etfit
     READ(iuna2Fsave,*) ((xkfit(i,ik), i=1,3), ik=1,nksfit)
     READ(iuna2Fsave,*) wkfit
     READ(iuna2Fsave,*) nk1fit, nk2fit, nk3fit
     CLOSE( UNIT = iuna2Fsave )
  END IF
  write(*,*) 'ionode', ionode    
  !!
  nkfit=nk1fit*nk2fit*nk3fit

  write (*,*) 'nksq, nksfit : ',nksq, nksfit

  allocate (eqkfit(nkfit), sfit(nkfit))
  allocate(vik(nkfit, nbnd, 3))
  !
  ! map k-points in the IBZ to k-points in the complete uniform grid
  !
  write(*,*) 'before lint fit grid, nkfit: ', nkfit
  nksfit_real=nksfit/nspin_lsda
  call lint ( nsym, s, .true., at, bg, npk, 0,0,0, &
       nk1fit,nk2fit,nk3fit, nksfit_real, xkfit, 1, nkfit, eqkfit, sfit)
  !! after this step, eqkfit is defined on a uniform fit grid with normal order
  !! eqkfit(a point on full fit grid) = a point in the reduced fit grid
  !! etfit is the energy defined for reduced fit grid points. eqkfit map full fit to reduced fit so can obtain energy from etfit
  deallocate (sfit, xkfit, wkfit)
  
 absdk(1) = sqrt(dot_product(bg(:,1)/nk1fit*2*3.14/(alat*0.53), bg(:,1)/nk1fit*2*3.14/(alat*0.53))) !! absolute |dk| in dir a, in ang-1
 absdk(2) = sqrt(dot_product(bg(:,2)/nk1fit*2*3.14/(alat*0.53), bg(:,2)/nk1fit*2*3.14/(alat*0.53)))
 absdk(3) = sqrt(dot_product(bg(:,3)/nk1fit*2*3.14/(alat*0.53), bg(:,3)/nk1fit*2*3.14/(alat*0.53)))

 vec(:,1) = bg(:,1)/nk1fit*2*3.14/(alat*0.53)/ absdk(1) !! unit vector in a
 vec(:,2) = bg(:,2)/nk1fit*2*3.14/(alat*0.53)/ absdk(2) !! unit vector in b
 vec(:,3) = bg(:,3)/nk1fit*2*3.14/(alat*0.53)/ absdk(3) !! unit vector in c
 
	ispin = 1
  !! compute velocity of full grid
  do ikz = 1, nk1fit
		do iky = 1, nk2fit
			do ikx = 1, nk3fit
				ik = (ikz -1) * nk1fit* nk2fit + (iky -1) * nk1fit + ikx
				ikxp1 =ikx +1
				ikyp1 =iky +1
				ikzp1 =ikz +1
				ikxm1 = ikx -1
				ikym1 = iky -1
				ikzm1 = ikz -1
				if(ikx .eq. nk1fit) then
					ikxp1 = 1
				endif
				if(iky .eq. nk2fit) then
					ikyp1 = 1
				endif
				if(ikz .eq. nk3fit) then
					ikzp1 = 1
				endif
				if(ikx .eq. 1) then
					ikxm1 = nk1fit
				endif
				if(iky .eq. 1) then
					ikym1 = nk2fit
				endif
				if(ikz .eq. 1) then
					ikzm1 = nk3fit
				endif
				indxp=(ikz -1) * nk1fit* nk2fit + (iky -1) * nk1fit + ikxp1 !! the three points to the right side of point ik
				indyp=(ikz -1) * nk1fit* nk2fit + (ikyp1 -1) * nk1fit + ikx !! used to compute velocity, finite difference, very brutal
				indzp=(ikzp1 -1) * nk1fit* nk2fit + (iky -1) * nk1fit + ikx
				indxm=(ikz -1) * nk1fit* nk2fit + (iky -1) * nk1fit + ikxm1 !! the three points to the right side of point ik
				indym=(ikz -1) * nk1fit* nk2fit + (ikym1 -1) * nk1fit + ikx !! used to compute velocity, finite difference, very brutal
				indzm=(ikzm1 -1) * nk1fit* nk2fit + (iky -1) * nk1fit + ikx
		  	do ibnd = bndsta, bndend
		  		etk = etfit(ibnd,eqkfit(ik)+nksfit*(ispin-1)/2)
		  		etknbp(1) = etfit(ibnd,eqkfit(indxp)+nksfit*(ispin-1)/2)
		  		etknbp(2) = etfit(ibnd,eqkfit(indyp)+nksfit*(ispin-1)/2)
		  		etknbp(3) = etfit(ibnd,eqkfit(indzp)+nksfit*(ispin-1)/2)
		  		etknbm(1) = etfit(ibnd,eqkfit(indxm)+nksfit*(ispin-1)/2)
		  		etknbm(2) = etfit(ibnd,eqkfit(indym)+nksfit*(ispin-1)/2)
		  		etknbm(3) = etfit(ibnd,eqkfit(indzm)+nksfit*(ispin-1)/2)
		  		do i=1, 3
		  			vik(ik, ibnd, i) = 0.0
		  			vtmp(:) = 0.5*(etknbp(:)-etknbm(:))*13.6/absdk(:)/(6.58e-6)  !! vik unit: m/s   1/hbar * (e+ - e-)/(2 dk) it's in the dir a, b, c
		  			do l = 1, 3
		  				vik(ik, ibnd, i) = vik(ik, ibnd, i) + vtmp(l) * vec(i, l)  !! now vik should be along x, y, z
		  			enddo
		  			if(abs(vik(ik, ibnd, i)) .lt. vthrs) then
		  				vik(ik, ibnd, i) =0.0d0 !! to correct for numerical error in computing v from finite difference, particularly when v small
		  			endif
! if(i .eq. 1) then
! 	write(145, '(2e14.4)')  etk*13.6, vik(ik, ibnd, i)
! endif
		  			!! central difference, not including / 2.0*dk, which doesn't matter		  			
		  		enddo
		  	enddo		  			  	
		  enddo
		enddo
	enddo


allocate (gf(nksq,nbnd,nsig, numnd), gf1(nksq,nbnd,nbnd, nsig, numnd))
gf = (0.0d0,0.0d0)
gf1 = (0.0d0,0.0d0)

!! read in wavefunctions
write(*,*) nksq
IF (nksq.GT.1) REWIND (unit = iunigk)
DO ik = 1, nksq
	IF (nksq.GT.1) THEN
		READ (iunigk, err = 100, iostat = ios) npw, igk
		100     CALL errore ('eimpmatrix', 'reading igk', ABS (ios) )
	ENDIF
	IF (lgamma) npwq = npw
	ikk = ikks(ik)  !! kpint ikk
	!! write(*,*) 'ikk, ik: ', ikk, ik
	IF (lsda) current_spin = isk (ikk)
	!
	! read unperturbed wavefuctions psi(k) and psi(k+q)
	!
	IF (nksq.GT.1) THEN
		IF (lgamma) THEN
			CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
		ELSE
			CALL errore ('eimpmatrix', 'lgamma not .true.!', ABS (ios) )
		ENDIF
	ENDIF
	do ibnd = bndsta, bndend !! only do a small range of bands
		CALL cft_wave (evc(1, ibnd), evr, +1)		
		evrtot(ik, ibnd-bndsta + 1, :, :) = evr(:,:)	!! FFT of all wavefunction into R space and store in evrtot
	enddo
ENDDO
write(*,*) 'Finished preparing wavfunctions'
!! read in energies
allocate(eig(nksq, nbnd))
filename = "./tmp/si.save/data-file.xml.eig"
open(93, file=filename, status="unknown")
do i= 1, 6 !! skip heading
	read(93, *) lines
enddo
cbm = 1.0d3
do ik = 1, nksq
	do i = 1, 4
		read(93, *) lines
	enddo
	do ibnd = 1, nbnd
		read(93, *) eig(ik, ibnd)
		eig(ik,ibnd) = eig(ik, ibnd)*2.0 !! in the unit of Ry, which is 2*Hartree
		if(ibnd == bndcbm .and. eig(ik,ibnd) .lt. cbm) then
			cbm = eig(ik,ibnd) !! finding cbm
			kcbm = ik
		endif
	enddo
	do i = 1, nbnd + 4
		read(93, *) lines !! skip occupations
	enddo
enddo
write(*,*) 'finished reading eigen energies'
close(93)

!! open(17, file="elimpmat.dat", status="unknown")	

xksft(:) = xk(:,kcbm)

write(*,*) 'CBM is found to be: ', cbm, ' Ry'
write(*,'(a12, i5, a1, 3f12.5)') 'the kpt is: ', kcbm, ' ', xksft(1), xksft(2), xksft(3)
if (parabo .eq. 1) write(*,*)  'PARABO:: shift kpoints to center at cbm'
if (parabo .eq. 1) cbm =0.0

thrs = thrs/13.6 !! convert from eV to Ry
!! compute rate
do ik=1,nksq			
	write(*,*) 'ik=', ik
	do ibnd = bndsta, bndend
		if(parabo .eq. 1) then
			ki(:)= (xk(:,ik)-xksft(:))*(2*3.14/alat)/bohr !! now in Ang
			kim = sqrt(dot_product(ki, ki)) !! length of kjm in Ang-1
			etk = pref*kim*kim !! Ry
		else 
			ki(:)= xk(:,ik)*(2*3.14/alat)/bohr !! now in Ang
			kim = sqrt(dot_product(ki, ki)) !! length of kjm in Ang-1
			etk = eig(ik, ibnd)
		endif
		if(abs(etk - cbm) .le. thrs) then	!! if the ik,ibnd is within the to-be-reported range
			v1(:) = vik(map(ik), ibnd, :)  !! velocity for ik, ibnd
			!! write(88, '(i8, 4e16.4)') ik, etk, v1(1), v1(2), v1(3) !! output the velocities around cbm for checking
			do jk = 1, nksq  !! do double counting						
				do jbnd = bndsta, bndend
					!! if (ik .eq. jk .and. ibnd .eq. jbnd ) cycle !! modes do not scatter into itself
					!! write(*,*) 'ik= ', ik, ' jk= ', jk				
					flag = 0
					if(parabo .eq. 1) then
						kj(:)= (xk(:,jk)-xksft(:))*(2*3.14/alat)/bohr !! now in Ang
						kjm = sqrt(dot_product(kj, kj)) !! length of kjm
						etkj = pref*kjm*kjm !! Ry
					else
						kj(:)= xk(:,jk)*(2*3.14/alat)/bohr !! now in Ang
						kjm = sqrt(dot_product(kj, kj)) !! length of kjm
						etkj = eig(jk, jbnd)
					endif
					if(abs(etkj - etk) .lt. 0.5/13.6) then  !! if state jk, jbnd is beyond 1eV from current state ik, ibnd, skip it			
						v2(:) = 		vik(map(jk), jbnd, :)	
						kimkj(:)= (xk(:,ik)- xk(:,jk))*(2*3.14/alat)/bohr !! cartesian coord of k_i - k_j , real coord: *2pi/alat	   ; 
						kimkjm=sqrt(dot_product(kimkj, kimkj)) !! units in 1/Ang
						do isig = 1, nsig
							!! write(*,*) 'ideg: ', isig, ', deg: ', deg(isig)
							degauss1 = deg(isig)
							!! thrswindow = 1.0/3.14* aimag(1.0/(3.0*degauss1 - cim*degauss1)) !! only beyond this energy will be counted as overlapping
							w0g1 = 1.0/3.14* aimag(1.0/(etk - etkj - cim*degauss1))							
							!! if (w0g1 .gt. thrswindow) then !! compute elimpmat only when needed
							if(flag .eq. 0) then
								if(dot_product(v1, v1) .lt. 1d-6 .or. dot_product(v2, v2) .lt. 1d-6) then
									costheta = -1 ! contribut
								else
									costheta = dot_product(v1,v2)/sqrt(dot_product(v1,v1)*dot_product(v2,v2))
								endif
								evir(:,:) = evrtot( ik, ibnd-bndsta + 1, :, :)
								evjr(:,:) = evrtot( jk, jbnd-bndsta + 1, :, :)															
								if (useakkp .eq. 1) then
									call compute_akkp(evir, evjr, dffts%nnr, 1,  akkp) !! compute Akk' =(1/V_cell)* integrate_cell(uk'* uk)
									!! write(49, *) abs(akkp)  !! check what Akk' I get
								else
									akkp = 1.0 !! test plane wave part
								endif
								do ind = 1, numnd
									call compute_elimpmat_analyt( Ld(ind), kimkjm, elimpplwv)
									!! write(60+ind,'(3e16.4)') kimkjm, real(elimpplwv), aimag(elimpplwv)  !! |<psi_i|dV|psi_j>|,  Ry*Ang^3
									eimpmt(ind) = akkp*elimpplwv !! Hkk' = Akk' * integrate(plwv_k Us plwv_k'); DOES NOT include (1/V_normalization) because it's absorbed outside Hkk'		
								enddo
								flag = 1 !! only need to compute once for different deg
							endif
							!! costheta = -1
							do ind = 1, numnd
								gf(ik, ibnd, isig, ind) = gf(ik, ibnd,isig, ind) + real(conjg(eimpmt(ind))* eimpmt(ind)) * wkinp *w0g1 *(1-costheta)			!! sum_k'jbnd |<k'j|dV|ki>|^2 * delta(e_k'j - e_ki)
								gf1(ik, ibnd, jbnd, isig, ind) = gf1(ik, ibnd, jbnd, isig, ind) + real(conjg(eimpmt(ind))* eimpmt(ind)) * wkinp *w0g1 *(1-costheta)
! 								if(ik .eq. kcbm .and. isig .eq. 4) then
! 									write(*,'(i5, 4e16.4, f6.2)') jk, w0g1, real(gf(ik, ibnd,isig, ind)), abs(eimpmt), wkinp, 1-costheta
! 								endif
							enddo
							!! endif
						enddo !! isig				
					endif
				enddo !! jk
			enddo ! ik
		endif
	enddo !! jbnd
enddo !! ibnd

do ik = 1, nksq
	do ibnd = bndsta, bndend
		if(parabo .eq. 1) then
			ki(:)= (xk(:,ik)-xksft(:))*(2*3.14/alat)/bohr !! now in Ang
			kim = sqrt(dot_product(ki, ki)) !! length of kjm in Ang-1
			etk = pref*kim*kim !! Ry
		else
			etk = eig(ik,ibnd)
		endif
		if( abs(etk - cbm) .le. thrs) then
			do ind = 1, numnd
				write(90+ind, '(2i6, 100e20.4)') ik, ibnd, etk*13.6, (real(gf(ik, ibnd, isig, ind))*2*3.14/(6.6d-16/13.6) *nd(ind)/((alat*bohr)**3/4.0*1.0d-30)* 1.0d-60, isig=1,nsig)    !! 0.25=1/4 is the primitive cell vol in cart, now units of s^-1
				!! if the interpolation of matrix elements is correct, this should be correct, only missing factor of n_d
				!! write(92, '(2i6, 2f20.4)') ik, ibnd, etk, nd* real(2.0*pi*gff(ik, ibnd, isig)) * onemcos * ryd2mev !!!!!!!!!!!!!!!! Only for debug fit=dense, so no interpolation is taking effect
				!! 1- cos_theta, theta the angle between two group velocities, here just take cos=0
				!! nd should relate to defect concentratio
			enddo
		endif
	enddo
enddo
do ind = 1, numnd
	do ik = 1, nksq
		etk = eig(ik,bndsta)
		if( abs(etk - cbm) .le. thrs) then
			write(100+ind, '(2i6, 100e20.4)') ik, bndsta, etk*13.6, (real(gf1(ik, bndsta,bndsta,  isig, ind))*2*3.14/(6.6d-16/13.6) *nd(ind)/((alat*bohr)**3/4.0*1.0d-30)* 1.0d-60, isig=1,nsig) 
			write(110+ind, '(2i6, 100e20.4)') ik, bndsta, etk*13.6, (real(gf1(ik, bndsta,bndend,  isig, ind))*2*3.14/(6.6d-16/13.6) *nd(ind)/((alat*bohr)**3/4.0*1.0d-30)* 1.0d-60, isig=1,nsig)
		endif		
		etkj = eig(ik,bndend)
		if( abs(etkj - cbm) .le. thrs) then
			write(120+ind, '(2i6, 100e20.4)') ik, bndend, etkj*13.6, (real(gf1(ik, bndend,bndsta,  isig, ind))*2*3.14/(6.6d-16/13.6) *nd(ind)/((alat*bohr)**3/4.0*1.0d-30)* 1.0d-60, isig=1,nsig)
			write(130+ind, '(2i6, 100e20.4)') ik, bndend, etkj*13.6, (real(gf1(ik, bndend,bndend,  isig, ind))*2*3.14/(6.6d-16/13.6) *nd(ind)/((alat*bohr)**3/4.0*1.0d-30)* 1.0d-60, isig=1,nsig)
		endif
	enddo
enddo

end subroutine compute_rate


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE compute_akkp(evir, evjr, nnr, npol, akkp)
use kinds, only : dp
implicit none
integer :: nnr, npol, ir
complex(dp) :: evir(nnr, npol), evjr(nnr, npol), akkp

akkp = 0.0d0
do ir = 1, nnr
	akkp = akkp + conjg(evjr(ir,1))*evir(ir,1)
enddo 
akkp = akkp / nnr

END SUBROUTINE compute_akkp

!-----------------------------------------------------------------------
SUBROUTINE compute_elimpmat_analyt( Ld, kimkjm, eimpmt)
!-----------------------------------------------------------------------
USE kinds, ONLY : DP
USE cell_base, ONLY : alat, at, bg
implicit none
real(DP) ::kimkjm, Ld
complex(dp), intent(out) :: eimpmt
real(dp) :: ks, e0, pi, mm, ee, hbar

	ks=11.68
	e0=8.854e-12
	pi=3.14
	mm=1.08*9.1e-31
	ee=1.6e-19
	hbar=1e-34

eimpmt = ee**2/(ks*e0)/(kimkjm*kimkjm + 1.0/Ld/Ld)/(1e-10)/(1.6e-19)/13.6
!! the prefix is ~ 1.13813809712

END SUBROUTINE compute_elimpmat_analyt