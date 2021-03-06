!
!        Copyright (C) 2000-2014 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM DS
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
subroutine COLLISIONS_sex(X,Xk,E,k,q,Xw)
 !
 ! SEX collisions 
 !
 use pars,          ONLY:SP,pi,cZERO
 use memory_m,      ONLY:mem_est
 use com,           ONLY:msg
 use frequency,     ONLY:w_samp
 use electrons,     ONLY:levels,spin_occ,spin
 use LIVE_t,        ONLY:live_timing
 use com,           ONLY:error
 use parallel_m,    ONLY:PP_redux_wait,PAR_IND_QP,PAR_IND_Q,&
&                        PAR_n_Bp_mat_elements,PAR_COM_QP_A2A,&
&                        PAR_IND_Bp_mat,PAR_nQP,PAR_nQ,PAR_COM_G_b_A2A,PAR_COM_Q_index
 use collision,     ONLY:elemental_collision,elemental_collision_free,&
&                        SEX_collisions,nSEX_collisions,COLL_bands,&
&                        COLLISIONS_IO_status,elemental_collision_alloc
 use IO_m,          ONLY:io_control,OP_RD_CL,REP,VERIFY,NONE,Io_and_Messaging_switch,io_COLLs
 use QP_m,          ONLY:QP_t,QP_n_G_bands,QP_nb,QP_n_states,QP_table
 use X_m,           ONLY:X_alloc,X_mat,X_t
 use interfaces,    ONLY:WF_load
 use wave_func,     ONLY:WF,WF_buffer,WF_buffered_IO
 use R_lattice,     ONLY:qindx_S,bz_samp
 use hamiltonian,   ONLY:B_mat_index
 !
 implicit none
 !
 type(levels)             ::E
 type(bz_samp)            ::k,q,Xk
 type(X_t)                ::X
 type(w_samp)             ::Xw
 !
 ! Work Space
 !
 integer                 :: i_qp,ig1,ig2,iqbz,iqibz,alloc_err,iqs,n_bands_to_load,&
&                           i_n,i_k,i_k_s,i_sp,i_m,i_mp,n_bands(2),N_g_vecs,&
&                           i_np,i_kmq,i_kmq_bz,i_kmq_s,i_kmq_G0,PAR_nQP_local
 !
 type(elemental_collision):: isc,iscp
 complex(SP)              :: pre_factor
 !
 complex(SP), allocatable :: EM1s(:,:)
 integer                  :: io_err,ID
 integer, external        :: io_X
 !
 ! Reset
 !
 alloc_err=0
 !
 call elemental_collision_free(isc)
 call elemental_collision_free(iscp)
 !
 ! Setup
 !
 call k_build_up_BZ_tables(k)
 call k_build_up_BZ_tables(q)
 call k_expand(k)
 !
 ! Dimensions
 !
 QP_n_G_bands=COLL_bands
 !
 n_bands   =QP_n_G_bands
 N_g_vecs  =X%ng
 !
 n_bands_to_load=max(QP_n_G_bands(2),QP_nb)
 !
 call msg('r', '[SEX] Plane waves :',N_g_vecs)
 !
 ! WFs
 !
 if (WF_buffered_IO) then
   !
   call WF_load(WF_buffer,0,1,COLL_bands,(/1,k%nibz/),space='B',title='-BUFFER')
   !
 endif
 !
 call WF_load(WF,N_g_vecs,maxval(qindx_S(:,:,2)),(/1,n_bands_to_load/),(/1,k%nibz/),title='-SEX Collisions')
 !
 isc%iqref =0
 !
 ! Plasmon-Pole/Static interaction DB I/O
 !
 call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1,2/),MODE=VERIFY,ID=ID)
 io_err=io_X(X,Xw,ID)
 if (io_err<0) call error('Incomplete and/or broken PPA/Static diel. fun. database')
 !
 ! Test the spatial Inversion
 !   
 call WF_spatial_inversion(E,Xk)
 !
 ! Allocation
 !
 call X_alloc('X',(/N_g_vecs,N_g_vecs,1/))
 allocate(EM1s(N_g_vecs,N_g_vecs),stat=alloc_err)
 call mem_est("EM1s",(/size(EM1s)/),errors=(/alloc_err/))
 !
 call elemental_collision_alloc(isc, NG=N_g_vecs,NG_GAMP=(/N_g_vecs,N_g_vecs/),TITLE="SEX")
 call elemental_collision_alloc(iscp,NG=N_g_vecs,TITLE="SEX")
 !
 ! Parallel Indexes & LIVE-TIMING STEPS
 !
 PAR_nQP_local=PAR_nQP
 if (allocated(COLLISIONS_IO_status)) then
   PAR_nQP_local=count(COLLISIONS_IO_status==-1)
 endif
! call live_timing("SEX collisions (e-e correlation)",PAR_nQP_local*PAR_nQ*PAR_n_Bp_mat_elements)
 call live_timing("SEX collisions (e-e correlation)",QP_n_states)
 !
 do i_qp=1,QP_n_states
   !
   if (.not.PAR_IND_QP%element_1D(i_qp)) cycle
   !
   if (io_COLLs) then
     if (COLLISIONS_IO_status(i_qp)/=-1) cycle
   endif
   !
   i_n     = QP_table(i_qp,1)
   i_m     = QP_table(i_qp,2)
   i_k     = QP_table(i_qp,3)
   i_sp    = spin(QP_table(i_qp,:))
   !
   SEX_collisions(1)%v3_c=cZERO
   !
   do iqbz=1,q%nbz
     !
     iqibz    =q%sstar(iqbz,1)
     iqs      =q%sstar(iqbz,2)
     !
     isc%qs   =(/0,iqibz,iqs/)
     iscp%qs  =(/0,iqibz,iqs/)
     !
     if (iqibz/=isc%iqref) then
       !
       call scatter_Gamp(isc,'c')
       !
       call Io_and_Messaging_switch("-fragments_IO_log")
       !
       ! I/O
       !
       if (q%nbz==1) then
         call io_control(ACTION=OP_RD_CL,COM=NONE, SEC=(/1,2,3/),ID=ID)
       else
         if (iqbz ==1) call io_control(ACTION=OP_RD_CL,COM=NONE, SEC=(/1,2,2*iqibz+1/),ID=ID)
         if (iqbz>  1) call io_control(ACTION=OP_RD_CL,COM=NONE, SEC=(/2*iqibz,2*iqibz+1/),ID=ID)
       endif
       io_err=io_X(X,Xw,ID)
       !
       call Io_and_Messaging_switch("+fragments_IO_log")
       !
       ! Poles and Residuals
       !
!$omp parallel do default(shared), private(ig1,ig2)
       do ig1=1,X%ng
         do ig2=1,X%ng
           EM1s(ig1,ig2)=X_mat(ig1,ig2,1)*isc%gamp(ig1,ig2)
         enddo
       enddo
!$omp end parallel do
       !
     endif
     !
     ! Scattering geometry
     !=====================
     !
     !                      W(i_kmq_G0,iqibz,iqs)                    
     !                        _ _ _ _ \ _ _ _ _
     !                       /        /        \
     !                      /                   \
     !                     /     G_{np,mp}       \
     !     (i_n,i_k,1)-->--.---------->----------.---->----(i_m,i_k,1)
     !                          (i_kmq,i_kmq_s)
     !
     i_k_s   = 1
     !
     i_kmq_bz =qindx_S(i_k,iqbz,1)
     i_kmq_G0 =qindx_S(i_k,iqbz,2)
     !
     i_kmq   = k%sstar(i_kmq_bz,1)
     i_kmq_s = k%sstar(i_kmq_bz,2)
     !
     isc%is  =(/i_n,i_k,  i_k_s,  i_sp/)
     isc%os  =(/0  ,i_kmq,i_kmq_s,i_sp/)
     !
     iscp%is =(/i_m,i_k,  i_k_s,  i_sp/)
     iscp%os =(/0,  i_kmq,i_kmq_s,i_sp/)
     !
     isc%qs(1) =i_kmq_G0
     iscp%qs(1)=i_kmq_G0
     !
     do i_mp=n_bands(1),n_bands(2)
       !
       iscp%os(1)=i_mp
       !
       call scatter_Bamp(iscp)
       !
       do i_np=n_bands(1),n_bands(2)
         !
         ! COLLISIONS(n,n',qp_n_states) = \sum_{G,G',q} ( rho(m,n,G,q,k) x EM1s(G ,G') x rho*(m',n',G',q,k) /|q+G|/|q+G'| )
         !
         isc%os(1)=i_np
         !
         if (any(isc%is/=iscp%is).or.any(isc%os/=iscp%os)) then
           call scatter_Bamp(isc)
         else
           isc%rhotw=iscp%rhotw
         endif
         !
         pre_factor=cZERO
         !
!$omp parallel do default(shared), private(ig1), reduction(+:pre_factor)
         do ig1=1,X%ng
           do ig2=1,X%ng
             pre_factor=pre_factor+isc%rhotw(ig2)*EM1s(ig2,ig1)*conjg(iscp%rhotw(ig1))
           enddo
         enddo
!$omp end parallel do
         !
         SEX_collisions(1)%v3_c(i_np,i_mp,i_kmq)=SEX_collisions(1)%v3_c(i_np,i_mp,i_kmq)-4._SP/spin_occ*pi*pre_factor
         !
       enddo ! loop on the band i_np
       !
     enddo ! loop on the band i_mp
     !
   enddo ! loop on iqbz
   !
   call live_timing(steps=1)
   !
   call COLLISIONS_linearize_and_IO("SEX",i_qp,SEX_collisions,nSEX_collisions)
   !
 enddo ! loop on i_qp
 !
 call live_timing()
 !
 ! CLEAN
 !
 deallocate(EM1s)
 call mem_est("EM1s")
 call X_alloc('X')
 call elemental_collision_free(isc)
 call elemental_collision_free(iscp)
 !
end subroutine
