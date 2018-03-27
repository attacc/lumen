! 
!        Copyright (C) 2017-2018 LUMEN team
!              http://www.attaccalite.com/lumen/ 
!
! Authors (see AUTHORS file for details): AC
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
module hamiltonians
 !
 use pars,       ONLY:SP,schlen,lchlen
 use electrons,  ONLY:levels
 !
 implicit none
 !
 integer           :: NE_bands(2)  ! Bands range
 integer           :: NE_nbands    ! Number of bands
 !
 integer           :: NE_ng      ! number of G-vectors in the Hartree/TDDFT
 integer           :: ng_oscll   ! number of G-vectors in the Oscillators
 ! 
 ! Initial Energies and Damping
 !------------
 real(SP),    allocatable :: E_full(:,:,:)
 !
 ! Bands in terms of Kohn-Sham orbitals
 complex(SP), allocatable :: V_bands(:,:,:,:) 
 !
 ! Time-dependent Valence bands in real-space
 complex(SP), allocatable :: VAL_BANDS(:,:,:,:)
 !
 ! Integrator Variables
 integer, parameter :: VER_HIGH = 1
 integer, parameter :: VER_LOW  = 2
 !
 integer :: NE_verbosity
 character(schlen)            :: Verbosity_name
 !
 ! Hamiltoniana flags
 ! 
 logical :: eval_HARTREE
 logical :: eval_DFT
 logical :: eval_OVERLAPS
 logical :: eval_dG
 !
 ! Response function Hamiltonian
 ! Subtract the equilibrium correlation part
 !
 logical :: l_resp_func_H 
 !
 ! Correlation flags
 !
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
end module
