!
!        Copyright (C) 2000-2016 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): CA AM DS
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
subroutine COLLISIONS_compose_nl(G_input)
 !
 ! HF(qp_state) = - (0.,1.) * \sum_{n,n',q} STATIC_collisions(n,n',q,qp_state) x G_lesser(n,n',k-q)
 ! 
 use pars,             ONLY:SP,cZERO,cI
 use drivers,          ONLY:l_sc_hf,l_sc_cohsex,l_sc_exx
 use QP_m,             ONLY:QP_table,QP_nk,QP_n_states
 use collision,        ONLY:Fock_collisions,SEX_collisions,nFock_collisions,nSEX_collisions
 use electrons,        ONLY:spin
 use parallel_m,       ONLY:PAR_QP_index,PAR_IND_QP,PP_redux_wait,PAR_COM_Xk_ibz_INDEX
 use nl_optics,        ONLY:NE_Vnl_xc,NL_bands
 !
 implicit none
 !
 complex(SP), intent(in)   :: G_input(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),QP_nk)
 !
 ! Work Space
 !
 logical           :: l_expand
 integer           :: i_qp,ib,ibp,i_k,i_kmq,ic1,ic2,i_qp_mem,i_spin
 integer           :: i_n,i_m
 !
 if(.not.l_sc_hf.and..not.l_sc_cohsex.and..not.l_sc_exx) return
 !
 NE_Vnl_xc=cZERO
 !
 do i_qp=1,QP_n_states
   !
   if (.not.(PAR_IND_QP%element_1D(i_qp))) cycle
   !
   i_qp_mem=PAR_QP_index(i_qp)
   !
   i_k=QP_table(i_qp,3)
   i_n=QP_table(i_qp,1)
   i_m=QP_table(i_qp,2)
   !
   ic1=1
   ic2=1
   !
   i_spin=spin(QP_table(i_qp,:))
   !
   do i_kmq=1,QP_nk
     !
     do ibp=NL_bands(1),NL_bands(2)
       do ib=NL_bands(1),NL_bands(2)
         !
         if (nFock_collisions>0) then
           if (Fock_collisions(i_qp_mem)%table(ib,ibp,i_kmq)/=0) then
             NE_Vnl_xc(i_n,i_m,i_k,i_spin)=NE_Vnl_xc(i_n,i_m,i_k,i_spin) &
&                 -cI*Fock_collisions(i_qp_mem)%v_c(ic1)*G_input(ib,ibp,i_kmq)
             ic1=ic1+1
           endif
         endif
         if (nSEX_collisions>0) then
           if (SEX_collisions(i_qp_mem)%table(ib,ibp,i_kmq)/=0) then
             NE_Vnl_xc(i_n,i_m,i_k,i_spin)=NE_Vnl_xc(i_n,i_m,i_k,i_spin) &
&                 -cI*SEX_collisions(i_qp_mem)%v_c(ic2)*G_input(ib,ibp,i_kmq)
             ic2=ic2+1
           endif
         endif
         !
       enddo
     enddo
     !
   enddo
   !
 enddo
 !
 ! Symmetrization
 !================
 do i_qp=1,QP_n_states
   !
   if (.not.(PAR_IND_QP%element_1D(i_qp))) cycle
   !
   i_k=QP_table(i_qp,3)
   i_n=QP_table(i_qp,1)
   i_m=QP_table(i_qp,2)
   !
   NE_Vnl_xc(i_n,i_n,i_k,:)=real(NE_Vnl_xc(i_n,i_n,i_k,:))
   NE_Vnl_xc(i_m,i_n,i_k,:)=conjg(NE_Vnl_xc(i_n,i_m,i_k,:))
   !
 enddo
 !
 call PP_redux_wait(NE_Vnl_xc,COMM=PAR_COM_Xk_ibz_INDEX%COMM)
 !
end subroutine COLLISIONS_compose_nl
