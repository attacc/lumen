!
!        Copyright (C) 2000-2017 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): CA
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
module lum_m
 !
 use drivers,    ONLY:l_lum
 use pars,       ONLY:SP
 !
 ! Luminescence Spectrum
 !
 real(SP),    allocatable :: LU(:,:)
 !
 ! Luminescence direction
 !
 real(SP)                 :: LU_weights(3)
 !
 ! BS dipoles
 !
 complex(SP), allocatable :: BSS_dipoles_LU(:,:)
 complex(SP), allocatable :: BSS_LU_f(:,:)
 !
 contains
   !
   subroutine LU_init( )
     !
     use pars,       ONLY:pi,rZERO
     use units,      ONLY:SPEED_OF_LIGHT
     use R_lattice,  ONLY:RL_vol,nqbz
     use vec_operate,ONLY:v_norm
     use BS_solvers, ONLY:BSS_n_freqs
     use memory_m,   ONLY:mem_est
     !
     implicit none
     !
     if (allocated(LU))  return
     !
     allocate(LU(BSS_n_freqs,3))
     LU = rZERO
     LU_weights=LU_weights/v_norm(LU_weights)
     !
   end subroutine
   !
   subroutine BSS_LU_alloc( )
     !
     use pars,       ONLY:cZERO,rZERO
     use BS,         ONLY:BS_H_dim,BS_K_dim
     use memory_m,   ONLY:mem_est
     !
     implicit none
     !
     if (.not.l_lum) return
     !
     allocate(BSS_dipoles_LU(BS_H_dim,3))
     call mem_est("BSS_dipoles_LU",(/size(BSS_dipoles_LU)/),(/2*SP/))
     BSS_dipoles_LU = cZERO
     !
     allocate(BSS_LU_f(2,BS_K_dim))
     call mem_est("BSS_LU_F",(/size(BSS_LU_F)/),(/SP/))
     BSS_LU_f = rZERO
     !
   end subroutine
   !
   subroutine BSS_LU_free( )
     !
     use memory_m,    ONLY:mem_est
     !
     implicit none
     !
     if (.not.l_lum) return
     !
     if (allocated(BSS_dipoles_LU)) then
       deallocate(BSS_dipoles_LU)
       call mem_est("BSE_dipoles_LU")
     endif
     if (allocated(BSS_LU_f)) then
       deallocate(BSS_LU_f)
       call mem_est("BSE_LU_F")
     endif
     !
   end subroutine
   !
end module lum_m 
