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
subroutine NL_Hamiltonian(E,k,q,X,Time,i_time,V_bands)
 !
 ! The subroutine calculate Ho + V + Sigma and the fields E_tot, E_ext, E_ks
 !
 !   H_k=Ho_k+U_k+ V_k^H-V_k^H[rho_ref]+V_xc[rho]
 !                                        
 use pars,           ONLY:SP,cZERO,pi,cONE
 use units,          ONLY:SPEED_OF_LIGHT
 use electrons,      ONLY:levels,n_sp_pol,n_spinor
 use R_lattice,      ONLY:bz_samp
 use zeros,          ONLY:zero_dfl
 use QP_m,           ONLY:QP_nk
 use X_m,            ONLY:X_t,global_gauge
 use fields,         ONLY:A_ext,A_vecpot,Efield
 use wave_func,      ONLY:WF
 use xc_functionals, ONLY:V_xc,XC_potential_driver
 use nl_optics,      ONLY:Correlation,V_xc_0,full_rho,IPA,E_full,LRC,JGM,LRCW,l_use_DIPOLES,NL_LRC_alpha,   &
&                         NL_initial_P,E_tot,E_ext,E_ks,E_xc_0,I_relax,Alpha_ED,use_OSCILLATORS,dG,eval_dG, & 
&                         V_hartree_ne,V_xc_ne,rho_n,Hk,eval_DFT,eval_HARTREE,rho_reference,NL_bands,NE_Vnl_xc,&
&                         Ho_plus_Sigma,NL_P,Phase_LifeTime,NL_nbands,VAL_BANDS,SEX,HF
 use global_XC,      ONLY:WF_KIND,WF_xc_functional 
 use parallel_m,     ONLY:PP_redux_wait,PAR_IND_Xk_ibz,PAR_COM_Xk_ibz_INDEX,master_cpu,myid,PAR_Xk_ibz_index
 use IO_m,           ONLY:IO_and_Messaging_switch
 use wrapper_omp,    ONLY:V_by_V_plus_V_omp
 !
 implicit none
 !
 type(levels),  intent(in)    :: E
 type(bz_samp), intent(in)    :: k,q
 real(SP),      intent(in)    :: Time
 type(X_t),     intent(in)    :: X
 complex(SP),   intent(in)    :: V_bands(NL_bands(2),E%nbf,QP_nk,n_sp_pol)
 integer,       intent(in)    :: i_time
 !
 ! Working Space
 !
 integer     :: ik,is,i1,ik_mem
 complex(SP) :: E_vec_pot(3),E_xc_pot(3)
 real(SP) :: alpha_
 !
 if(eval_HARTREE .or.  eval_DFT) then
   !      
   call NL_build_valence_bands(E,V_bands,VAL_BANDS,NL_bands(1)) 
   !
   call el_density_vbands(E,k,rho_n,VAL_BANDS)
   !
 endif
 !
 if(eval_HARTREE) call V_Hartree(rho_n-rho_reference,V_hartree_ne)
 !
 ! I reconstruct the full density if necessary because both TDDFT 
 ! and JGM are not linear in rho_n and require the full density
 ! 
 if(NL_bands(1)>1.and.Correlation/=IPA) rho_n=full_rho+rho_n-rho_reference
 !
 if(eval_DFT) then
   !
   ! Reconstruct the full density
   !
   call XC_potential_driver(E,k,WF_KIND,WF_xc_functional,1,MUTE="-report -log -io_in -io_out",unMUTE=master_cpu)
   call IO_and_Messaging_switch("+report +log +io_in +io_out")
   V_xc_ne=V_xc-V_xc_0
   !
 elseif(Correlation/=IPA) then
   V_xc_ne=cZERO
 endif
 !
 ! Evaluate the Polarization 
 ! ===================================================================
 !
 if(l_use_DIPOLES.or.trim(global_gauge)=='velocity') then
   call DIP_polarization(E,k,V_bands,NL_P)
 else
   ! This subroutine update also the Sm1_plus,Sm1_minus and All_S_det matrices
   call Berry_Polarization_NEQ(E,k,V_bands,NL_P)
   !
 endif
 !
 NL_P=NL_P-NL_initial_P ! Subtract the initial polarization
 !
 ! Evaluate the external/total field
 !=================================
 !
 call NL_ext_fields(A_ext,Time)
 !
 A_vecpot=A_ext
 !
 E_ext=-A_ext%vecpot_vel/SPEED_OF_LIGHT
 !
 E_tot=E_ext
 E_ks =E_ext
 !
 if(Correlation==LRC.or.Correlation==JGM.or.Correlation==LRCW) then
    !
    ! Evaluate the induced field and/or the LRC correction
    !
    E_vec_pot=-A_ext%vecpot_vel/SPEED_OF_LIGHT
    !
    E_tot=E_vec_pot
    !
    if(Correlation==LRCW) then !This can be calculated once for each frequency, but I did not want to mess up NL_driver (MG)
      alpha_ = NL_LRC_alpha 
      if (Efield(1)%ef_name=='SOFTSIN'.or.Efield(1)%ef_name=='SIN') &
      &     alpha_  = Alpha_ED(E,Efield(1)%frequency)
           E_vec_pot=E_vec_pot + alpha_*NL_P
    endif
    !
    if(Correlation==LRC)   E_vec_pot=E_vec_pot + NL_LRC_alpha*NL_P
    if(Correlation==JGM) then
      call EXC_macroscopic_JGM(E,NL_P,rho_n,E_xc_pot,V_xc_ne)
      E_vec_pot=E_vec_pot +  E_xc_pot - E_xc_0
    end if
    !
    E_ks=E_vec_pot
    !
    A_vecpot%vecpot_vel=-E_vec_pot*SPEED_OF_LIGHT
    !
 endif
 ! 
 if(eval_dG)   call NL_build_dG_lesser(E,V_bands,dG)
 ! 
 ! Evaluate collisions if requires
 !=======================
 if((Correlation==SEX.or.Correlation==HF).and..not.use_OSCILLATORS) then
   call COLLISIONS_compose_nl(dG)
 endif
 !
 !
 ! Build the Hamiltonian
 ! =====================
 !
 Ho_plus_Sigma=cZERO
 I_relax      =cZERO
 !
 do ik=1,QP_nk
   !
   if (.not.PAR_IND_Xk_ibz%element_1D(ik)) cycle
   ik_mem=PAR_Xk_ibz_index(ik)
   !
   ! Here I fill Hk with...
   !
   Hk=cZERO
   !
   ! Correlation part Hartree + TDDFT
   ! ================================
   !
   if(Correlation/=IPA) then
     !  
     if(n_spinor==1) then
        call V_real_space_to_H(ik,V_hartree_ne+V_xc_ne,Hk,WF,'def')
     elseif (n_spinor==2) then
        call V_real_space_to_H(ik,V_hartree_ne,Hk,WF,'def')
        if(eval_DFT) call V_real_space_to_H(ik,V_xc_ne,Hk,WF,'xc ')
     endif
     !
   endif
   !
   if((Correlation==SEX.or.Correlation==HF).and.use_OSCILLATORS) then
     call NL_SEX(V_bands,Hk,k,q,E,ik)
   endif
   !
   ! Add Ho = T+V_ion+V_h+V_xc
   !==================
!$omp parallel do default(shared), private(i1,is)
   do i1=NL_bands(1),NL_bands(2)
     Hk(i1,i1,:)=Hk(i1,i1,:)+E_full(i1,ik_mem,:)
   enddo
!$omp end parallel do
   ! 
   if(l_use_DIPOLES.or.trim(global_gauge)=='velocity') then
     !    
     call NL_apply_field(ik,Hk,A_vecpot)
     !
   else
     !
     ! Use Sm1_plus, Sm1_minus and V_bands to build the W operator
     !
     do is=1,n_sp_pol
       call  Build_W_operator(V_bands,E,k,A_vecpot,ik,is,Hk)
     enddo
     !
   endif                                                         
   !
   if((Correlation==SEX.or.Correlation==HF).and..not.use_OSCILLATORS) then
     call V_by_V_plus_V_omp(NL_nbands,cONE,NE_Vnl_xc(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),ik,1),&
&                                                 Hk(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),1))

   endif
   !
   ! Damping: NL_damping fills the I_relax matrix 
   ! ================================================
   call NL_damping(E,V_bands,ik)
   !
   ! Store the Hamiltonian
   ! =====================
!$omp parallel do default(shared), private(i1)
   do i1=NL_bands(1),NL_bands(2)
     Ho_plus_Sigma(NL_bands(1):NL_bands(2),i1,ik_mem,:)=Hk(NL_bands(1):NL_bands(2),i1,:)
   enddo
!$omp end parallel do
   !
 enddo  !---- loop on k-point
 !
end subroutine NL_Hamiltonian

