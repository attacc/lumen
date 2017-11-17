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
subroutine NL_driver(E,X,k,q,Xw)
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
 use QP_m,           ONLY:QP_nk
 use QP_CTL_m,       ONLY:QP_apply
 use frequency,      ONLY:w_samp
 use nl_optics,      ONLY:n_runs,Correlation,IPA,NL_alloc,NL_free,V_bands,NL_bands,NL_step,&
&                         Integrator,l_use_DIPOLES,NL_P_t,NL_ng,EF_parms,Sigma_SEX_EQ,ng_oscll,   &
&                         E_ext_t,E_tot_t,E_ks_t,E_ks,E_tot,E_ext,runs_done,n_runs_todo,LSEX,LHF, &
&                         NL_nbands,NL_P,NL_i_time,NL_steps,NL_alloc_real_space_quantities, &
&                         NL_alloc_k_distributed,E_bare
 use IO_m,           ONLY:IO_and_Messaging_switch
 use interfaces,     ONLY:PARALLEL_global_indexes,PARALLEL_WF_distribute,PARALLEL_WF_index
 use parallel_m,     ONLY:PAR_IND_runs,PAR_IND_runs_ID,CPU_structure,PAR_IND_Xk_ibz,PAR_COM_Xk_ibz_INDEX, &
&                         PP_redux_wait,myid,PAR_Xk_ibz_index 
 !
 implicit none
 !
 ! Notice that E is never changed!
 !
 type(levels)      :: E
 type(bz_samp)     :: k,q
 type(X_t)         :: X(2)
 type(w_samp)      :: Xw(2)
 logical, external :: stop_now
 ! 
 ! Work Space
 !
 integer           :: i_r,i1,NL_bands1_save,N_runs_save,iGo_max,ik,ik_mem
 real(SP)          :: Time
 integer, external :: eval_G_minus_G
 integer           :: wf_ng_save
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
 ! QP dimensions
 !================
 call QP_setup(E)
 !
 ! Parallelization
 ! =======================
 call setup_DIPOLES_parallelization() 
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
 call check_non_periodic()
 !
 NL_bands(1)=NL_bands1_save
 !
 ! Free memory: unused DIPOLES
 ! ===========================
 !
 call X_alloc('P_square')
 if(l_use_DIPOLES) then
   call X_alloc('DIP_S')
 else
   call X_alloc('DIP_R')
   call X_alloc('DIP_P')
 endif
 !
 ! Allocation
 ! ==================
 call NL_alloc(E,q) 
 !
 ! Evaluate Oscillators for LSEX/LHF
 ! =============================
 if(Correlation==LSEX.or.Correlation==LHF) then
   !
   ! Parallelization on k and q
   !
   call setup_OSCL_parallelization()
   !
   ! Evaluate oscillators and write them on disk
   call NL_oscillators_eval(q,k,E)
   !
   X(2)%ng=ng_oscll
   call NL_lsex_potential(q,X(2),Xw(2)) 
   !
 else
   !      
   ! Free memory: qindx_* 
   ! ====================
   call qindx_free()
   !
 endif
 !
 ! Generate DIPOLES in all the BZ if required
 !===================
 if(l_use_DIPOLES) call Build_DIP_bz(k,E,X)
 !
 ! Non_Linear parallization
 ! ===================================
 call setup_NL_parallelization() 
 !
 ! Reload WF with less G-vecs for the real-time dynamics
 ! =====================================================
 !
 wf_ng_save=wf_ng
 wf_ng     =NL_ng
 iGo_max   =1
 call WF_free(WF)
 call PARALLEL_WF_distribute(K_index=PAR_IND_Xk_ibz,CLEAN_UP=.TRUE.)
 !
 ! WF loading
 !=============
 if(Correlation/=IPA) then
   !  
   call WF_load(WF,NL_ng,iGo_max,(/1,NL_bands(2)/),(/1,k%nibz/),space='R',title='-NL')  
   !
   ! Allocate all quantities that depends from the fft_size
   ! ======================================================
   call NL_alloc_real_space_quantities(E,q)
   !
 endif
 !
 ! Allocation with k-distributed memory
 ! ==================
 call NL_alloc_k_distributed(E) 
 !
 ! QP-corrections 
 !================
 call QP_apply(NL_bands,E,k,"G",msg_fmt='s')
 !
 ! Load oscillators
 ! ===================
 if(Correlation==LSEX.or.Correlation==LHF) call NL_oscillators_load(q,k)
 !    
 ! Initialization
 ! ==================
 call NL_start_and_restart(E,k,X)
 !
 ! Equilibrium Exchange
 !=========================
 ! If LSEX or LHF calculate the equilibrium contribution 
 ! to the screened exchange
 !
 if(Correlation==LSEX.or.Correlation==LHF) then
   Sigma_SEX_EQ=cZERO
   do ik=1,QP_nk
     if (.not.PAR_IND_Xk_ibz%element_1D(ik)) cycle
     ik_mem=PAR_Xk_ibz_index(ik)
     call NL_long_range_SEX(V_bands,Sigma_SEX_EQ(:,:,:,ik_mem),k,q,E,ik)
   enddo
 endif
 !
 ! Check if calculation was already done
 ! =====================================
 if(n_runs_todo==0) return
 !
 ! Initialize IO
 ! ==================
 call NL_databases_IO(0,'write')
 !
 ! Start-Up & Timing
 !===================
 call live_timing('[NL] Runs :',PAR_IND_runs%n_of_elements(PAR_IND_runs_ID+1)*NL_steps)
 !
 do i_r=1,n_runs
   !
   if(.not.PAR_IND_runs%element_1D(i_r)) cycle
   !
   if(runs_done(i_r)) then
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
   Efield(1)%frequency=EF_parms(i_r)%frequency
   Efield(1)%intensity=EF_parms(i_r)%intensity
   Efield(1)%amplitude=EF_parms(i_r)%amplitude
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
     E_ext_t(NL_i_time,:)=real(E_ext)
     E_tot_t(NL_i_time,:)=real(E_tot)
     E_ks_t (NL_i_time,:)=real(E_ks)
     !
     call NL_output(NL_i_time,i_r) 
     !
     call NL_Integrator(E,k,q,X,Integrator,NL_i_time,V_bands)
     !
     call live_timing(steps=1)
     !
   enddo
   !
   call NL_databases_IO(i_r,'write')
   !
 enddo
 !
 call live_timing()
 !
 ! Deallocation
 ! ==========================
 !
 wf_ng=wf_ng_save
 !
 call WF_free(WF)
 call NL_free()
 call X_alloc('DIP_S')
 !
 call set_real_printed_length()
 !
 contains
   !
   subroutine setup_DIPOLES_parallelization()
     use parallel_m,       ONLY:CPU_and_ROLE_strings_save
     use nl_optics,        ONLY:eval_OVERLAPS
     implicit none
     !
     ! Set bands between 1 and NL_bands(2)
     !
     NL_bands1_save  =NL_bands(1)
     NL_bands(1)     =1
     X(1)%ib         =(/1,NL_bands(2)/)
     !
     ! Set n_runs = 1 
     ! Parallelization on frequencies turned off
     !
     N_runs_save=n_runs
     n_runs   =1
     eval_OVERLAPS   =.TRUE.
     !
     ! Save user-defined parallel structure
     ! and force parallelization on k-points (for dipoles and overlaps)
     call CPU_and_ROLE_strings_save()
     !
     call PARALLEL_global_indexes(E,k,q," ",RESET=.TRUE.,X=X(1))
     call PARALLEL_global_indexes(E,k,q,"Non_Linear",RESET=.FALSE.,X=X(1))
     !
   end subroutine setup_DIPOLES_parallelization
   !
   subroutine setup_OSCL_parallelization
     use parallel_m,       ONLY:PAR_IND_Xk_ibz,CPU_and_ROLE_strings_restore, &
 &                              PAR_IND_k_plus_q
     use nl_optics,        ONLY:eval_OVERLAPS
     implicit none
     !
     ! Restore original number of bands
     !
     X(1)%ib       =(/NL_bands(1),NL_bands(2)/)
     !
     ! Restore user-defined parallel structure
     call CPU_and_ROLE_strings_restore()
     !
     call PARALLEL_global_indexes(E,k,q," ",RESET=.TRUE.,X=X(1))
     call PARALLEL_global_indexes(E,k,q,"Oscillators",X=X(1))
     call PARALLEL_WF_distribute(K_index=PAR_IND_k_plus_q,CLEAN_UP=.TRUE.)
     call PARALLEL_WF_index( )
     !
   end subroutine setup_OSCL_parallelization
   !
   subroutine setup_NL_parallelization
     use parallel_m,       ONLY:PAR_IND_Xk_ibz,CPU_and_ROLE_strings_restore
     use nl_optics,        ONLY:eval_OVERLAPS
     implicit none
     !
     ! Restore original number of bands
     !
     X(1)%ib       =(/NL_bands(1),NL_bands(2)/)
     !
     ! Restore frequencies
     ! Parallelization on frequencies turned on
     !
     n_runs =N_runs_save
     eval_OVERLAPS   =.FALSE.
     !
     ! Restore user-defined parallel structure
     call CPU_and_ROLE_strings_restore()
     !
     call PARALLEL_global_indexes(E,k,q," ",RESET=.TRUE.,X=X(1))
     call PARALLEL_global_indexes(E,k,q,"Non_Linear",X=X(1))
     call PARALLEL_WF_distribute(K_index=PAR_IND_Xk_ibz,CLEAN_UP=.TRUE.)
     call PARALLEL_WF_index( )
     !
   end subroutine setup_NL_parallelization
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
   subroutine check_non_periodic()
    use R_lattice,      ONLY:k_map
    use D_lattice,      ONLY:a
    use com,            ONLY:error
    use fields,         ONLY:Efield
    use zeros,          ONLY:zero_norm
    implicit none
    !
    integer :: id_red
    !
    if(all(k_map%max_kdir/=1)) return
    !
    do id_red=1,3
      if(k_map%max_kdir(id_red)/=1) cycle
      if(abs(sum(a(id_red,:)*Efield(1)%versor(:)))>zero_norm) &
          call error(" Non-periodic directions not implemented yet ! ")
    enddo
    !
   end subroutine check_non_periodic
   !
end subroutine NL_driver
