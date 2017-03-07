!
!        Copyright (C) 2000-2015 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): MG CA
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
subroutine NL_driver(E,X,k,q)
 !
 use pars,           ONLY:SP,cZERO,cONE
 use stderr,         ONLY:set_real_printed_length
 use com,            ONLY:error
 use LIVE_t,         ONLY:live_timing 
 use electrons,      ONLY:levels
 use X_m,            ONLY:X_t,Dipole_approach,X_alloc
 use wave_func,      ONLY:WF,wf_ng
 use interfaces,     ONLY:WF_load
 use R_lattice,      ONLY:bz_samp,qindx_free
 use memory_m,       ONLY:mem_est
 use fields,         ONLY:Efield
 use QP_CTL_m,       ONLY:QP_apply
 use R_lattice,      ONLY:qindx_S
 use nl_optics,      ONLY:n_frequencies,Correlation,IPA,NL_alloc,NL_free,V_bands,NL_bands, &
&                         NL_estep,NL_er,Integrator,l_use_DIPOLES,NL_P_t,  &
&                         E_ext_t,E_tot_t,E_ks_t,E_ks,E_tot,E_ext,freqs_done,n_freqs_todo, &
&                         NL_nbands,NL_P,NL_i_time,NL_steps,NL_step
 use IO_m,           ONLY:IO_and_Messaging_switch
 use interfaces,     ONLY:PARALLEL_global_indexes,PARALLEL_WF_distribute,PARALLEL_WF_index
 use parallel_m,     ONLY:PAR_IND_WF_k,PAR_IND_freqs,PAR_IND_freqs_ID,CPU_structure
 !
 implicit none
 !
 ! Notice that E is never changed!
 !
 type(levels)      :: E
 type(bz_samp)     :: k,q
 type(X_t)         :: X(2)
 logical, external :: stop_now
 ! 
 ! Work Space
 !
 integer           :: i_f,i1,NL_bands1_save,N_frequency_save
 real(SP)          :: Time
 integer, external :: eval_G_minus_G
 !
 ! The non-linear bands
 ! ==============================================
 !
 NL_nbands=NL_bands(2)-NL_bands(1)+1
 !
 call section('+','Initialization')
 !#################################
 !
 call set_real_printed_length(f_length=10,g_length=10)  
 !
 ! Initialize
 !=========
 call NL_initialize(E,k,q)
 !
 ! Deallocate qindx_S (used only in previous subroutine)
 ! ========================
 if(allocated(qindx_S)) then
   deallocate(qindx_S)
   call mem_est("qindx_S")
 endif
 !
 ! QP dimensions
 !================
 call QP_setup(E)
 !
 ! Parallelization
 ! =======================
 NL_bands1_save  =NL_bands(1)
 NL_bands(1)     =1
 X(1)%ib         =(/1,NL_bands(2)/)
 N_frequency_save=n_frequencies
 n_frequencies   =1
 call PARALLEL_global_indexes(E,k,q," ",RESET=.TRUE.,X=X(1))
 call PARALLEL_global_indexes(E,k,q,"Non_Linear",RESET=.FALSE.,X=X(1))
 !
 ! Dipoles
 !=========
 !
 ! before loading the WF in R-space I check and/or calculate
 ! the optical collisions that are done in G-space and would
 ! force the reloading of the WFs
 !
 Dipole_approach='Covariant' 
 if(l_use_DIPOLES) then
   call DIPOLE_driver(E,k,X(1),X(1)%q0)
 else
   if(load_Overlaps()/=0) call DIPOLE_driver(E,k,X(1),X(1)%q0)
 endif
 !
 NL_bands(1)=NL_bands1_save
 !
 ! Free memory: qindx_* and unused DIPOLES
 ! =====================================================
 call qindx_free()
 !
 call X_alloc('P_square')
 if(l_use_DIPOLES) then
   call X_alloc('DIP_S')
 else
   call X_alloc('DIP_R')
   call X_alloc('DIP_P')
 endif
 !
 ! Non_Linear parallization
 ! ===================================
 X(1)%ib       =(/NL_bands(1),NL_bands(2)/)
 n_frequencies =N_frequency_save
 call PARALLEL_global_indexes(E,k,q," ",RESET=.TRUE.,X=X(1))
 call PARALLEL_global_indexes(E,k,q,"Non_Linear",RESET=.FALSE.,X=X(1))
 call PARALLEL_WF_distribute(K_index=PAR_IND_WF_k,CLEAN_UP=.TRUE.)
 call PARALLEL_WF_index( )
 !
 ! WFs loading
 !=============
 if(Correlation/=IPA) &
&  call WF_load(WF,wf_ng,1,(/1,NL_bands(2)/),(/1,k%nibz/),space='R',title='-NL')  
 !
 ! QP-corrections 
 !================
 call QP_apply(NL_bands,E,k,"G",msg_fmt='s')
 !
 ! Allocation
 ! ==================
 call NL_alloc(E) 
 !
 ! Generate DIPOLES in all the BZ if required
 !===================
 if(l_use_DIPOLES) call Build_DIP_bz(k,E,X)
 !
 !
 ! Initialization
 ! ==================
 call NL_start_and_restart(E,k,X)
 !
 ! Check if calculation was already done
 ! =====================================
 if(n_freqs_todo==0) return
 !
 ! Initialize IO
 ! ==================
 call NL_databases_IO(0,'write')
 !
 ! Start-Up & Timing
 !===================
 call live_timing('[NL] Frequencies :',PAR_IND_freqs%n_of_elements(PAR_IND_freqs_ID+1)*NL_steps)
 !
 do i_f=1,n_frequencies
   !
   if(.not.PAR_IND_freqs%element_1D(i_f)) cycle
   !
   if(freqs_done(i_f)) then
     call live_timing(steps=NL_steps)
     cycle
   endif
   !
   ! Re-initialization
   ! ================
   V_bands=cZERO
   forall(i1=1:E%nbf) V_bands(i1,i1,:,:)=cONE
   !
   ! Efield
   ! ================
   Efield(1)%frequency=NL_er(1)+real(i_f-1,SP)*NL_estep
   !
   do NL_i_time=1,NL_steps
     !
     Time=real(NL_i_time-1,SP)*NL_step
     !
     call NL_Hamiltonian(E,k,q,X,Time,NL_i_time,V_bands)
     !
     ! Store Variables
     !
     NL_P_t (NL_i_time,:)=real(NL_P,SP)
     E_ext_t(NL_i_time,:)=E_ext
     E_tot_t(NL_i_time,:)=E_tot
     E_ks_t (NL_i_time,:)=E_ks
     !
     call NL_output(NL_i_time,i_f) 
     !
     call NL_Integrator(E,k,q,X,Integrator,NL_i_time,V_bands)
     !
     call live_timing(steps=1)
     !
   enddo
   !
   call NL_databases_IO(i_f,'write')
   !
 enddo
 !
 call live_timing()
 !
 ! Deallocation
 ! ==========================
 !
 call WF_free(WF)
 call NL_free()
 call X_alloc('DIP_S')
 !
 call set_real_printed_length()
 !
 contains
   !
   integer function load_Overlaps()
     !
     use pars,          ONLY:SP
     use X_m,            ONLY:DIP_S
     use IO_m,           ONLY:io_control,OP_RD_CL,VERIFY,REP
     use electrons,      ONLY:n_sp_pol
     !
     implicit none
     !
     ! Word Space
     !
     integer :: ID_S,io_err
     integer, external :: io_Overlaps
     integer :: ik,i_sp_pol
     !
     call X_alloc('DIP_S',(/X(1)%ib(2),X(1)%ib(2),6,k%nbz/))
     DIP_S=cZERO
     !
     io_err=0
     !
     call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1/),MODE=VERIFY,ID=ID_S)
     io_err=io_Overlaps(X,ID_S)
     !
     if(io_err==0) then
       call live_timing('[NL] Overlaps loading :',k%nbz)
       do ik=1,k%nbz
         do i_sp_pol=1,n_sp_pol
           call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1+ik+(i_sp_pol-1)*k%nbz/),ID=ID_S)
           io_err=io_err+io_Overlaps(X,ID_S)
         enddo
         call live_timing(steps=1)
       enddo
       call live_timing()
     endif
     !
     load_Overlaps=io_err
     !
   end function load_Overlaps
   !
   subroutine QP_setup(E)
     use R_lattice,      ONLY:nkibz
     use electrons,      ONLY:levels
     use QP_m,           ONLY:QP_nk,QP_state,QP_nb
     use memory_m,       ONLY:mem_est
     !
     implicit none
     !
     type(levels) :: E
     !
     ! Build-up QP table
     !
     QP_nb=NL_bands(2)
     QP_nk=nkibz
     ! 
     ! QP_state my be read from DBs at start-up
     !
     if (allocated(QP_state)) deallocate(QP_state)
     !
     allocate(QP_state(QP_nb,QP_nk))
     call mem_est("QP_state",(/QP_nb*QP_nk/))
     QP_state=.TRUE.
     call QP_state_table_setup(E)
     !
   end subroutine QP_setup
   !
end subroutine NL_driver
