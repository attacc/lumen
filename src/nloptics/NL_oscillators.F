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
subroutine NL_oscillators(q,k,en) 
 !
 ! Evaluate oscillators for the Local-SEX
 !
 use nl_optics,           ONLY:NL_bands,OSCL,W,ng_oscll
 use LIVE_t,              ONLY:live_timing
 use pars,                ONLY:SP,cZERO
 use R_lattice,           ONLY:bz_samp,qindx_S
 use electrons,           ONLY:levels
 use wave_func,           ONLY:WF
 use interfaces,          ONLY:WF_load
 use QP_m,                ONLY:QP_nk
 use collision,           ONLY:elemental_collision,elemental_collision_free, &
&                              elemental_collision_alloc         
 use parallel_m,          ONLY:PAR_IND_Xk_ibz,PAR_IND_Xk_ibz_ID
 use FFT_m,               ONLY:fft_size
 !
 implicit none
 !
 type(bz_samp), intent(in) :: q,k
 type(levels),  intent(in) :: en
 !
 ! Work Space
 !
 type(elemental_collision)            :: isc 
 integer :: i_k,i_m,i_n,i_q,iqbz,i_kmq_G0,iqs,i_kmq,i_kmq_s,i_k_s
 integer :: i_k_sp,i_kmq_sp,iqibz,ig,i_np,i_mp,io,ic
 !
 call elemental_collision_free(isc)
 !
 isc%iqref =0
 isc%ngrho =ng_oscll  ! only G=1 is required for local-SEX
 !
 call elemental_collision_alloc(isc,NG=isc%ngrho,TITLE="OSC")
 !
 call live_timing('[NL] Oscillators :',PAR_IND_Xk_ibz%n_of_elements(PAR_IND_Xk_ibz_ID+1))
 !
 do i_k=1,QP_nk
   do iqbz=1,q%nbz
   !
   i_kmq_G0  =qindx_S(i_k,iqbz,2)
   iqibz     =q%sstar(iqbz,1)
   iqs       =q%sstar(iqbz,2)
   isc%qs  =(/i_kmq_G0,iqibz,iqs/)
   !
   i_k_s=1
   i_kmq   = k%sstar(qindx_S(i_k,iqbz,1),1)
   i_kmq_s = k%sstar(qindx_S(i_k,iqbz,1),2)
   !
   !  Spin-Polarization not supported!!
   ! 
   i_k_sp  = 1
   i_kmq_sp= 1
   !
   do i_n=NL_bands(1),NL_bands(2)
     do i_m=NL_bands(1),NL_bands(2)
       isc%is  =(/i_n,i_k,i_k_s,i_k_sp/)
       isc%os  =(/i_m,i_kmq,i_kmq_s,i_kmq_sp/)
       call scatter_Bamp(isc)            
       OSCL(i_n,i_m,i_k,iqbz,:)=isc%rhotw(:)
       write(100,*) abs(OSCL(i_n,i_m,i_k,iqbz,1))
     enddo
   enddo
   !
   enddo
   !
   call live_timing(steps=1)
   !
 enddo
 !
 call live_timing()
 !
end subroutine NL_oscillators
