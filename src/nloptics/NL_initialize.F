!
!        Copyright (C) 2000-2015 the YAMBO team
!            http://www.yambo-code.org
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
subroutine NL_initialize(E,k,q)
 ! 
 use pars,           ONLY:SP,DP
 use drivers,        ONLY:l_sc_hf,l_sc_cohsex
 use vec_operate,    ONLY:normalize_v
 use units,          ONLY:AU2VMm1,AUT2FS,HA2EV
 use electrons,      ONLY:levels,n_sp_pol
 use com,            ONLY:msg,warning,error
 use memory_m,       ONLY:mem_est
 use zeros,          ONLY:zero_dfl
 use R_lattice,      ONLY:bz_samp,bare_qpg
 use fields,         ONLY:Efield,Efield_strength,n_ext_fields,EtoT
 use SC,             ONLY:SC_potential
 use collision,      ONLY:COLL_bands
 !
 ! PDE integrators
 ! ===============
 use nl_optics,      ONLY:EULER,EULEREXP,RK2,RK4,RK2EXP,RK4EXP,HEUN,INVINT,Integrator,CRANK_NICOLSON,Integrator_name
 !
 use nl_optics,      ONLY:NL_damping,NL_correlation,Correlation,NL_er,NL_estep,n_frequencies,IPA,TDDFT,HARTREE,LRC,JGM, &
&                         LRCW,eval_DFT,eval_HARTREE,NL_bands,NL_nbands,l_use_DIPOLES,SEX,HF,n_intensities, &
&                         NL_steps,NL_tot_time,Phase_LifeTime,NL_step,Verbosity_name,VER_HIGH,VER_LOW,NL_verbosity, &
&                         use_COLLISIONS,eval_dG
 !
 use electric,       ONLY:l_force_SndOrd
 use parser_m,       ONLY:parser
 !
 implicit none
 !
 type(levels), intent(in) :: E
 type(bz_samp),intent(in) :: q,k
 !
 ! Working space
 !
 integer  :: i1
 !
 ! Parser
 ! ============
 !
 call parser('UseDipoles',           l_use_DIPOLES)
 call parser('FrSndOrd',             l_force_SndOrd)
 !
 if(l_use_DIPOLES) call warning(' Fixed dipoles: only linear response properties are correct!')
 !
 ! Check double-precision
 ! ========================
 !
 if(DP/=SP)      call error(" Non-linear optics works only in double precision, configure with --enable-dp flag")
 if(n_sp_pol/=1) call error(" Spin-polarized systems not full implemented")
 if(any(NL_bands(:)==0)) call error(" Error set NLBands /= 0 ! ")
 if(NL_bands(1)>E%nbf.or.NL_bands(2)<=E%nbf) call error(" Error in bands range! ")
 !
 ! Integrator ! ============
 !
 select case (trim(Integrator_name))
    case('EULER')
        Integrator=EULER
    case('EULEREXP')
        Integrator=EULEREXP
    case('RK2')  
        Integrator=RK2
    case('RK4')
        Integrator=RK4
    case('RK2EXP')
        Integrator=RK2EXP
    case('RK4EXP')
        Integrator=RK4EXP
    case('HEUN')
        Integrator=HEUN
    case('INVINT')
        Integrator=INVINT
    case('CRANKNIC')
        Integrator=CRANK_NICOLSON
    case default
        call error("Unknow integrator")
 end select
 !
 ! Verbosity Level
 ! 
 select case (trim(Verbosity_name))
    case('low')
        NL_verbosity=VER_LOW
    case('high')
        NL_verbosity=VER_HIGH
    case default
        call error("Unknow verbosity level")
 end select
 !
 ! Important in the Time-dependent Schrodinger Eq.
 ! the Phase_LifeTime is 1/2 of the one used in the 
 ! Non-Equilibrium Green Functions
 !
 Phase_LifeTime=2._SP/NL_damping
 !
 ! Determine the total time from the damping
 !
 if(NL_tot_time<=0._SP) then
   !
   ! Run lenght is 6*Phase_lifeTime + sampling time
   !
   if(NL_damping<zero_dfl) call error(" Error damping=0 and NLtime<0!! ") 
   !
   ! Dephasing Time
   !
   NL_tot_time=Phase_LifeTime*5._SP
   !
   ! Add time to sample the response
   !
   NL_tot_time=NL_tot_time+EtoT(E=NL_er(1))+5._SP*NL_step
   !
 else
   call msg('s',' Total simulation time read from input') 
 endif
 !
 call msg('sr','Dephasing Time          [fs]:',5._SP*Phase_LifeTime*AUT2FS)
 call msg('sr','Sampling  Time          [fs]:',(EtoT(E=NL_er(1))+5._SP*NL_step)*AUT2FS)
 call msg('sr','Total simulation time   [fs]:',NL_tot_time*AUT2FS)
 !
 NL_steps=int(NL_tot_time/NL_step)+1
 !
 ! Correlation
 ! ========================
 !
 eval_HARTREE   =.false.
 eval_DFT       =.false.
 eval_dG        =.false.
 use_COLLISIONS =.false.
 !
 select case (trim(NL_Correlation))
    case('IPA')
        Correlation=IPA
        SC_potential='IP'
    case('HARTREE')
        Correlation=HARTREE
        eval_HARTREE =.true.
        SC_potential='HARTREE'
    case('TDDFT')  
        Correlation=TDDFT
        eval_HARTREE =.true.
        eval_DFT     =.true.
        SC_potential='default'
    case('LRC')
        Correlation=LRC
        SC_potential='HARTREE'
        eval_HARTREE =.true.
    case('LRW')
        Correlation=LRCW
        SC_potential='HARTREE'        
        eval_HARTREE =.true.
    case('JGM')
        Correlation=JGM
        SC_potential='HARTREE'        
        eval_HARTREE =.true.
    case('SEX')
        Correlation=SEX
        SC_potential   ='COHSEX'        
        eval_dG      =.true.
        eval_HARTREE =.true.
        l_sc_cohsex=.true.
    case('HF')
        Correlation=HF
        SC_potential   ='FOCK'
        eval_dG      =.true.
        eval_HARTREE =.true.
        l_sc_hf=.true.
    case default
        call error("Unknow correlation in the response function")
 end select
 !
 call parser('UseCollisions',use_COLLISIONS)
 !
 if(Correlation/=HF.and.Correlation/=SEX) use_COLLISIONS=.false.
 !
 !Quasi-particle setup
 !====================
 !
 if(use_COLLISIONS) then
   COLL_bands=NL_bands
   call QP_state_table_setup(E)
 endif
 !
 ! Deallocate bare_qpg if not required
 !=============================
 if(.not.eval_HARTREE.and.CORRELATION/=JGM) then
   deallocate(bare_qpg)
   call mem_est("bare_qpg")
 endif
 !
 ! Frequencies range
 ! =========================
 NL_estep=(NL_er(2)-NL_er(1))/real(n_frequencies,SP)
 !
 ! Generate Electric fields list
 ! ==============================
 call NL_generate_ef_list(Efield(1))
 !
 ! Intensities range
 !=========================
 if(n_intensities>1) &
&   call msg( 'sr', 'Multiple intensities for Richardson extrapolation')
 !
 ! Check Field
 ! =========================
 !
 ! This part of the code works with only one field
 !
 do i1=1,n_ext_fields
   !
   Efield(i1)%ef_name=trim(Efield(i1)%ef_name)      
   !
   if((Efield(i1)%ef_name=='none')) cycle
   !
   if(Efield(i1)%intensity==0._SP) &
&    call error(" External "//trim(Efield(i1)%ef_name)//" field has zero intensity.")
     !
   if(all(abs(Efield(i1)%versor(:))<epsilon(1.))) &
&    call error(" External "//trim(Efield(i1)%ef_name)//" field has undefined direction.")
     !
     Efield(i1)%versor        =normalize_v(Efield(i1)%versor)
     Efield(i1)%amplitude     =Efield_strength(Efield(i1)%intensity,'AU')
     if(Efield(i1)%t_initial_indx<2) then
       call warning(' Starting index should be >= 2 !')
     endif             
     Efield(i1)%t_initial     =real(Efield(i1)%t_initial_indx-1)*NL_step
     Efield(i1)%ef_pol        ="linear"
     !
     call check_symmetries(Efield(i1)%versor)
     !
 enddo
 !
 if(Efield(1)%ef_name/='SOFTSIN'.and.Efield(1)%ef_name/='SIN') &
    call warning(' Non-linear optics properties can be obtained with SOFTSIN/SIN fields only!')
 !
 ! LOG
 ! ===
 !
 call section('=','Messaging')
 !
 call msg( 'r','---General---')
 call msg( 'r', 'Integrator                 :',trim(Integrator_name))
 call msg( 'r', 'Simulation time        [fs]:',real(NL_steps-1)*NL_step*AUT2FS)
 call msg( 'r', 'Step length            [fs]:',NL_step*AUT2FS)
 call msg( 'r', 'NE Steps                   :',NL_steps)
 call msg( 'r', 'Correlation                :'//trim(NL_correlation))
 call msg( 'r', 'Damping                [eV]:',NL_damping*HA2EV)
 call msg( 'r', 'Dephasing              [fs]:',Phase_LifeTime*AUT2FS)
 call msg( 'r', 'Efield                     :'//trim(Efield(1)%ef_name))
 call msg( 'r', 'Efield damping         [fs]:',Efield(1)%damping*AUT2FS)
 call msg( 'r', 'Efield Amplitude      [V/m]:',Efield(1)%amplitude*AU2VMm1)
 call msg( 'r', 'Efield Initial Time    [fs]:',Efield(1)%t_initial*AUT2FS)
 call msg( 'r', 'Frequencies range          :',NL_er*HA2EV)
 call msg( 'sr','Frequencies steps          :',n_frequencies)
 call msg( 'sr','Intensity steps            :',n_intensities)
 call msg( 'r', 'Frequency   step           :',NL_estep*HA2EV)
 call msg( 'r', 'Use Dipoles                :',l_use_DIPOLES)
 call msg( 'r', 'Use COLLISIONS             :',use_COLLISIONS)
 call msg( 'r', 'Force Snd Order Dipoles    :',l_force_SndOrd)
 !
end subroutine NL_initialize
