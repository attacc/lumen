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
subroutine NL_build_valence_bands(en,V_bands,Val_bands) 
 !
 ! Construct the density matrix differences dGamma = gamma(t) - gamma(t=0) 
 ! from the time-dependent valence bands
 !
 use pars,        ONLY:SP,cZERO
 use electrons,   ONLY:levels,spin_occ,n_sp_pol
 use nl_optics,   ONLY:NL_bands
 use FFT_m,       ONLY:fft_size
 use wrapper_omp, ONLY:Vstar_dot_V_omp
 use QP_m,        ONLY:QP_nk
 use parallel_m,  ONLY:PP_redux_wait,PAR_IND_Xk_ibz,PAR_COM_Xk_ibz_INDEX
 implicit none
 !
 type(levels), intent(in)  :: en
 complex(SP),  intent(in)  :: V_bands(NL_bands(2),en%nbf,QP_nk,n_sp_pol)
 complex(SP),  intent(out) :: Val_bands(fft_size,en%nbf,QP_nk,n_sp_pol)
 !
 ! Work Space
 !
 integer :: i1,i2
 integer :: ik
 !
 Val_bands=cZERO
 !
 do ik=1,QP_nk
   !
   if (.not.PAR_IND_Xk_ibz%element_1D(ik)) cycle
   !
 enddo
 !
 call PP_redux_wait(Val_bands,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
 !
end subroutine NL_build_valence_bands
