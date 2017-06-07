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
 use nl_optics,           ONLY:NL_bands,OSCLL,W,ng_oscll
 use LIVE_t,              ONLY:live_timing
 use pars,                ONLY:SP,cZERO
 use R_lattice,           ONLY:bz_samp,qindx_S
 use electrons,           ONLY:levels
 use wave_func,           ONLY:WF,wf_ng
 use interfaces,          ONLY:WF_load
 use QP_m,                ONLY:QP_nk
 use collision,           ONLY:elemental_collision,elemental_collision_free, &
&                              elemental_collision_alloc         
 use parallel_m,          ONLY:PAR_IND_Xk_ibz,PAR_IND_Xk_ibz_ID
 use FFT_m,               ONLY:fft_size
 use IO_m,                ONLY:io_control,OP_RD_CL,VERIFY,REP,OP_WR_CL,OP_APP_WR_CL 
 use com,                 ONLY:msg
 use parallel_m,          ONLY:PP_redux_wait,PAR_IND_Xk_ibz,PAR_COM_Xk_ibz_INDEX,master_cpu
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
 integer :: i_k_sp,i_kmq_sp,iqibz,ig,i_np,i_mp,io,ic,i_sp_pol
 integer :: ID,io_err,iGo_max
 integer, external :: io_OSCLL
 !
 OSCLL=cZERO
 !
 i_sp_pol=1
 call msg('s','[OSC] Checking oscillators header ')
 call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1/),MODE=VERIFY,ID=ID)
 io_err=io_OSCLL(q,ID)
 !
 call PP_redux_wait(COMM=PAR_COM_Xk_ibz_INDEX%COMM ) 
 ! 
 if(io_err==0) then
   !
   ! Load Oscillators
   !   
   do i_k=1,QP_nk
     call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1+i_k+(i_sp_pol-1)*k%nibz/),ID=ID)
     io_err=io_OSCLL(q,ID)
   enddo
   !
   call msg('rs','[OSC] Oscillators loaded')
   !
   return
   !
 else
   !
   ! Oscillators not found, recalculated them
   !
   if(master_cpu) then
     call io_control(ACTION=OP_WR_CL,COM=REP,SEC=(/1/),MODE=VERIFY,ID=ID)
     io_err=io_OSCLL(q,ID)
   endif
   !
 endif
 !
 iGo_max=maxval(qindx_S(:,:,2)) 
 call WF_load(WF,wf_ng,iGo_max,(/1,NL_bands(2)/),(/1,k%nibz/),space='R',title='-NL')  
 !
 call elemental_collision_free(isc)
 !
 isc%iqref =0
 isc%ngrho =ng_oscll  ! only G=1 is required for local-SEX
 !
 call elemental_collision_alloc(isc,NG=isc%ngrho,TITLE="OSCLL")
 !
 call live_timing('[OSC] Oscillators :',PAR_IND_Xk_ibz%n_of_elements(PAR_IND_Xk_ibz_ID+1))
 !
 do i_k=1,QP_nk
   !
   if (.not.PAR_IND_Xk_ibz%element_1D(i_k)) cycle
   ! 
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
       OSCLL(i_n,i_m,i_k,iqbz,:)=isc%rhotw(:)
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
 do ig=1,ng_oscll
   call PP_redux_wait(OSCLL(:,:,:,:,ig),COMM=PAR_COM_Xk_ibz_INDEX%COMM ) 
 enddo
 !
 ! Write Oscillators on disk
 !
 do i_k=1,QP_nk
   !
   i_sp_pol=1
   call io_control(ACTION=OP_APP_WR_CL,COM=REP,SEC=(/1+i_k+(i_sp_pol-1)*k%nibz/),ID=ID)
   io_err=io_OSCLL(q,ID)
   !
 enddo
 !
end subroutine NL_oscillators
