! 
!        Copyright (C) 2000-2015 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): MG AC
! 
! This file is distributed under the terms of the GNU 
! General Public License. You can redistribute it and/or 
! modify it under the terms of the GNU General Public 
! License as published by the Free Software Foundation; 
! either version 2, or (at your option) any later version.
!
! This program is distributed in the hope that it will 
! be useful, but WITHOUT ANY WARRANTY; without even the 
! implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License 
! for more details.
!
! You should have received a copy of the GNU General Public 
! License along with this program; if not, write to the Free 
! Software Foundation, Inc., 59 Temple Place - Suite 330,Boston, 
! MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
!
module nl_optics
 !
 use pars,       ONLY:SP,schlen,lchlen
 use electrons,  ONLY:levels
 !
 implicit none
 !
 integer           :: NL_bands(2)  ! Bands range
 integer           :: NL_nbands    ! Number of bands
 !
 ! Long range correction field
 !
 real(SP)          :: NL_LRC_alpha
 !
 !
 integer           :: NL_steps  ! Number of real-time steps
 integer           :: NL_i_time ! Integer time-step
 real(SP)          :: NL_step   ! time-step in atomic units of time
 real(SP)          :: NL_tot_time ! Total simulation time
 !
 !
 integer           :: NL_ng      ! number of G-vectors in the Hartree/TDDFT
 integer           :: ng_oscll   ! number of G-vectors in the Oscillators
 ! 
 ! Polarization 
 ! ------------
 ! Calculated from the string-averaged discretized geometric phase
 ! see ref. PRB 47, 1651 (1993), or from DIPOLES
 !
 complex(SP)    :: NL_initial_P(3) ! in principle it can be different from zero
 !
 complex(SP)    :: NL_P(3)
 !
 ! Current
 ! 
 complex(SP)    :: NL_J(3)
 !
 ! Diagonal (QP) and off-diagonal (Phase) relaxation time
 !
 real(SP)    :: NL_damping
 real(SP)    :: Phase_LifeTime ! dephasing time
 !
 ! Initial Energies and Damping
 !------------
 complex(SP), allocatable :: I_relax(:,:,:,:,:)
 real(SP),    allocatable :: E_full(:,:,:)
 !
 ! For TDDFT, V_xc_0 stores the V_xc at equilibrium
 !
 real(SP), allocatable    :: V_xc_0(:,:)
 !
 ! For JGM E_xc_0 for the microscopic part of the functional
 !
 complex(SP) :: E_xc_0(3)
 !
 ! WFs and Density Matrix
 !
 ! Bands in terms of Kohn-Sham orbitals
 complex(SP), allocatable :: V_bands(:,:,:,:) 
 !
 ! Time-dependent Valence bands in real-space
 complex(SP), allocatable :: VAL_BANDS(:,:,:,:)
 !
 ! Integrator Variables
 !
 integer, parameter :: EULER    = 1
 integer, parameter :: EULEREXP = 2
 integer, parameter :: RK2      = 3
 integer, parameter :: RK4      = 4
 integer, parameter :: RK2EXP   = 5
 integer, parameter :: RK4EXP   = 6
 integer, parameter :: HEUN     = 7
 integer, parameter :: INVINT   = 8
 integer, parameter :: CRANK_NICOLSON = 9
 !
 integer                      :: Integrator
 character(schlen)            :: Integrator_name
 !
 ! Verbosity Level
 !
 integer, parameter :: VER_HIGH = 1
 integer, parameter :: VER_LOW  = 2
 !
 integer :: NL_verbosity
 character(schlen)            :: Verbosity_name
 !
 ! Correlation
 !
 character(schlen)            :: NL_correlation
 integer                      :: Correlation
 !
 integer, parameter :: IPA    = 1
 integer, parameter :: HARTREE= 2
 integer, parameter :: TDDFT  = 3
 integer, parameter :: LRC    = 4
 integer, parameter :: JGM    = 5
 integer, parameter :: LRCW   = 6
 integer, parameter :: SEX   = 7
 integer, parameter :: HF    = 8
 !
 logical :: eval_HARTREE
 logical :: eval_DFT
 logical :: eval_OVERLAPS
 logical :: eval_dG
 !
 ! Non-linear spectra parameters
 !
 real(SP) :: NL_er(2)
 real(SP) :: NL_estep
 integer  :: n_frequencies
 integer  :: n_intensities
 integer  :: n_runs_todo
 real(SP) :: initial_Amplitude
 integer  :: n_runs
 !
 integer, parameter :: n_parms=2 ! Number of paramters in efield loop
 !
 type EF_parameters
   real(SP) :: frequency
   real(SP) :: intensity
   real(SP) :: amplitude
   integer  :: indexes(n_parms) 
 end type
 !
 type(EF_parameters), allocatable :: EF_parms(:)
 !
 ! Databases input/output
 !
 integer                  :: ID(3)
 !
 ! Other
 !
 logical  :: l_use_DIPOLES   ! Use standard dipoles (valid only for linear response)
 !
 ! IO variables
 !
 real(SP),     allocatable :: NL_P_t(:,:)
 real(SP),     allocatable :: E_tot_t(:,:),E_ext_t(:,:), E_ks_t(:,:)
 complex(SP)               :: E_tot(3),E_ext(3),E_ks(3)
 logical,      allocatable :: runs_done(:)
 integer                   :: n_NL_descriptions
 character(lchlen)         :: NL_description(100)
 !
 ! Hamiltonian and potentials
 !
 complex(SP), allocatable :: Ho_plus_Sigma(:,:,:,:)   ! Full Hamiltonian
 real(SP),    allocatable :: E_bare(:,:,:)            ! Bare energies
 complex(SP), allocatable :: Hk(:,:,:)
 real(SP),    allocatable :: V_xc_ne(:,:)
 real(SP),    allocatable :: V_hartree_ne(:,:)
 !
 ! Oscillators
 !
 logical :: use_OSCILLATORS
 logical :: l_test_oscll
 complex(SP), allocatable :: OSCLL(:,:,:,:,:) ! Oscillators
 complex(SP), allocatable :: OSCLL_k(:,:,:,:) ! Oscillators in k, just for the evaluation
 complex(SP), allocatable :: W(:,:,:)
 complex(SP), allocatable :: Sigma_SEX_EQ(:,:,:,:)
 !
 ! Collisions
 !
 logical :: use_COLLISIONS
 complex(SP), allocatable :: dG(:,:,:) 
 complex(SP), allocatable :: NE_Vnl_xc(:,:,:,:)
 !
 ! Density
 !
 real(SP),    allocatable :: rho_reference(:)    ! reference density from all bands  
 real(SP),    allocatable :: full_rho(:)       ! density of all bands from 1:NL_bands(2)
 real(SP),    allocatable :: rho_n(:)          ! density of all bands from 1:NL_bands(2)
 !
 contains
   !
   subroutine NL_alloc_real_space_quantities(en,q)
     use FFT_m,          ONLY:fft_size
     use pars,           ONLY:SP
     use memory_m,       ONLY:mem_est
     use electrons,      ONLY:levels,n_spin,n_sp_pol
     use R_lattice,      ONLY:nXkbz,bz_samp
     use QP_m,           ONLY:QP_nk
     use xc_functionals, ONLY:V_xc
     !
     implicit none
     !
     type(levels),  intent(in) :: en
     type(bz_samp), intent(in) :: q
     !
     ! Allocate all quantities that depend from the fft_size
     !
     allocate(rho_reference(fft_size))
     call mem_est("rho_reference",(/fft_size/),(/SP/))
     !
     allocate(full_rho(fft_size))
     call mem_est("full_rho",(/fft_size/),(/SP/))
     !
     allocate(rho_n(fft_size))
     call mem_est("rho_n",(/fft_size/),(/SP/))
     !
     if(Correlation==TDDFT.or.Correlation==JGM) then
       allocate(V_xc_0(fft_size,n_spin))
       call mem_est("V_xc_0",(/fft_size,n_spin/),(/SP/))
       allocate(V_xc(fft_size,n_spin))
       call mem_est("V_xc",(/fft_size,n_spin/),(/SP/))
     endif          
     !
     allocate(V_xc_ne(fft_size,n_spin))
     call mem_est('V_xc_ne',(/size(V_xc_ne)/),(/SP/))
     !
     allocate(V_hartree_ne(fft_size,n_spin))
     call mem_est("Hartree_Potential",(/size(V_hartree_ne)/))
     !
   end subroutine NL_alloc_real_space_quantities
   !
   subroutine NL_alloc(en,q)
     !
     use pars,           ONLY:SP
     use memory_m,       ONLY:mem_est
     use electrons,      ONLY:levels,n_spin,n_sp_pol
     use R_lattice,      ONLY:nXkbz,bz_samp
     use QP_m,           ONLY:QP_nk
     use electric,       ONLY:ELECTRIC_alloc
     !
     implicit none
     !
     type(levels),  intent(in) :: en
     type(bz_samp), intent(in) :: q
     !
     allocate(Hk(NL_bands(2),NL_bands(2),n_sp_pol))
     call mem_est("Hk",(/size(Hk)/),(/2*SP/))
     !
     allocate(V_bands(1:NL_bands(2),en%nbf,QP_nk,n_sp_pol))
     call mem_est("V_bands",(/size(V_bands)/),(/2*SP/))
     !
     if(Correlation==SEX.or.Correlation==HF) then
       allocate(W(q%nbz,ng_oscll,ng_oscll))
       call mem_est("W",(/size(W)/),(/2*SP/))
     endif          
     !
     if(.not.l_use_DIPOLES) then
       call  ELECTRIC_alloc(en)
     endif
     !
     allocate(NL_P_t(NL_steps,3))
     call mem_est("NL_P_t",(/size(NL_P_t)/))
     !
     allocate(E_tot_t(NL_steps,3),E_ext_t(NL_steps,3),E_ks_t(NL_steps,3))
     call mem_est("E_tot_t E_ext_t",(/size(E_tot_t),size(E_ext_t),size(E_ks_t)/))
     !
     if((Correlation==SEX.or.Correlation==HF).and.use_COLLISIONS) then
       !  
       allocate(dG(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),QP_nk))
       call mem_est("dG",(/size(dG)/),(/2*SP/))
       !
       allocate(NE_Vnl_xc(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),QP_nk,n_sp_pol))
       call mem_est("NE_Vnl_xc",(/size(NE_Vnl_xc)/),(/2*SP/))
       NE_Vnl_xc=(0._SP,0._SP)
       !
     endif
     !
   end subroutine NL_alloc
     !
   subroutine NL_alloc_k_distributed(en)
     !  
     use parallel_m,     ONLY:PAR_Xk_nibz,PAR_Xk_ibz_index,PAR_nQP
     use pars,           ONLY:SP
     use memory_m,       ONLY:mem_est
     use electrons,      ONLY:levels,n_spin,n_sp_pol
     use R_lattice,      ONLY:nXkbz
     use QP_m,           ONLY:QP_nk
     use FFT_m,          ONLY:fft_size
     !
     implicit none
     !
     type(levels),  intent(in) :: en
     !
     ! Valence bands in fft space
     !
     allocate(VAL_BANDS(fft_size,en%nbf,PAR_Xk_nibz,n_sp_pol))
     call mem_est("Val_bands",(/size(VAL_BANDS)/),(/2*SP/))
     !
     if((Correlation==SEX.or.Correlation==HF).and.use_OSCILLATORS) then
       allocate(Sigma_SEX_EQ(NL_bands(2),NL_bands(2),n_sp_pol,PAR_Xk_nibz))
       call mem_est("Sigma_SEX_EQ",(/size(Sigma_SEX_EQ)/),(/2*SP/))
     endif          
     !
     allocate(I_relax(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),en%nbf,PAR_Xk_nibz,n_sp_pol))
     call mem_est("I_relax",(/size(I_relax)/),(/2*SP/))
     !
     allocate(Ho_plus_Sigma(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),PAR_Xk_nibz,n_sp_pol))
     call mem_est("Ho_plus_Sigma",(/size(Ho_plus_Sigma)/),(/2*SP/))
     !
     ! Reference energies
     !
     allocate(E_bare(1:NL_bands(2),PAR_Xk_nibz,n_sp_pol))
     call mem_est("E_bare",(/size(E_bare)/))
     E_bare=en%E(:NL_bands(2),:,:)
     !
     allocate(E_full(1:NL_bands(2),PAR_Xk_nibz,n_sp_pol))
     call mem_est("E_full",(/size(E_full)/))
     !
   end subroutine NL_alloc_k_distributed
   
   subroutine NL_free
     !
     use memory_m,       ONLY:mem_est
     use electric,       ONLY:ELECTRIC_free
     use xc_functionals, ONLY:V_xc
     !
     implicit none
     !
     deallocate(Ho_plus_Sigma)
     call mem_est("Ho_plus_Sigma")
     !
     deallocate(Hk)
     call mem_est("Hk")
     !
     deallocate(I_relax)
     call mem_est("I_relax")
     !
     deallocate(V_bands)
     call mem_est("V_bands")
     !
     if(allocated(VAL_BANDS)) then
       deallocate(VAL_BANDS)
       call mem_est("Val_bands")
     endif
     !
     if(allocated(rho_n)) then
       deallocate(rho_n)
       call mem_est("rho_n")
     endif
     !
     if(allocated(rho_reference)) then
       deallocate(rho_reference)
       call mem_est("rho_reference")
     endif
     !
     if(allocated(full_rho)) then
       deallocate(full_rho)
       call mem_est("full_rho")
     endif
     !
     if(allocated(E_bare)) then
         deallocate(E_bare)
         call mem_est("E_bare")
     endif
     !
     if(allocated(E_full)) then
       deallocate(E_full)
       call mem_est("E_full")
     endif
     !
     if(allocated(V_xc_0)) then
       deallocate(V_xc_0)
       call mem_est("V_xc_0")
     endif
     !
     if(allocated(V_xc)) then
       deallocate(V_xc)
       call mem_est("V_xc")
     endif
     !
     if(allocated(V_xc_ne)) then
       deallocate(V_xc_ne)
       call mem_est("V_xc_ne")
     endif
     !
     if(allocated(V_hartree_ne)) then
       deallocate(V_hartree_ne)
       call mem_est("V_hartree_ne")
     endif
     !
     if(.not.l_use_DIPOLES) then
       call ELECTRIC_free()
     endif
     !
     deallocate(NL_P_t)
     call mem_est("NL_P_t")
     !
     deallocate(E_tot_t,E_ext_t,E_ks_t)
     call mem_est("E_tot_t E_ext_t E_ks_t")
     !
     if(Correlation==SEX.or.Correlation==HF) then
       if(use_OSCILLATORS) then
         deallocate(OSCLL)
         call mem_est("OSCLL")
         deallocate(Sigma_SEX_EQ)
         call mem_est("Sigma_SEX_EQ")
         deallocate(W)
         call mem_est("W")
       endif
       if(use_COLLISIONS) then
         deallocate(dG)
         call mem_est("dG")
         deallocate(NE_Vnl_xc)
         call mem_est("NE_Vnl_xc")
       endif
     endif
     !
   end subroutine NL_free
   !
   real(SP) function Alpha_ED(en,omega)
   !
   real(SP) :: Egap,beta,omega
   type(levels) :: en
     !	
     Egap = sum(en%E(en%nbf+1,:,1)-en%E(en%nbf,:,1))/real(en%nk,SP)
     beta = NL_LRC_alpha/Egap**2
     Alpha_ED = NL_LRC_alpha + beta*omega**2
   end function
   !
end module
