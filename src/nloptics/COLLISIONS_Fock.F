!
!        Copyright (C) 2000-2015 the YAMBO team
!              http://www.yambo-code.org
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
subroutine COLLISIONS_Fock(E,k,q,i_k,Fock_Collisions)
 !
 ! Fock Collisions
 !
 use pars,          ONLY:SP,schlen,pi,cZERO
 use electrons,     ONLY:levels,spin_occ,spin
 use collision,     ONLY:elemental_collision,elemental_collision_free,elemental_collision_alloc
 use QP_m,          ONLY:QP_ng_Sx
 use R_lattice,     ONLY:qindx_S,bz_samp
 use nl_optics,     ONLY:NL_bands
#if defined _TIMING
 use timing_m,      ONLY:timing
#endif
 use wave_func,     ONLY:WF 
 use interfaces,    ONLY:WF_load
 use wrapper_omp,   ONLY:V_dot_V_omp
 !
 implicit none
 !
 type(levels)             ::E
 type(bz_samp)            ::k,q
 !
 integer,  intent(in)     :: i_k
 complex(SP), intent(out) :: Fock_collisions(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),k%nibz)
 !
 ! Work Space
 !
 integer                 :: i_qp,ig1,iqbz,iqibz,alloc_err,iqs,&
&                           i_n,i_k_s,i_sp,i_m,i_mp,n_bands(2),N_g_vecs,&
&                           i_np,i_kmq,i_kmq_bz,i_kmq_s,i_kmq_G0,PAR_nQP_local
 !
 type(elemental_collision)           :: isc,iscp
 complex(SP)             :: pre_factor
 !
 ! Reset
 !
 call elemental_collision_free(isc)
 call elemental_collision_free(iscp)
 !
 ! Dimensions
 !
 n_bands   =NL_bands
 N_g_vecs  =QP_ng_Sx
 !
 ! WFs
 !
 call WF_load(WF,N_g_vecs,maxval(qindx_S(:,:,2)),n_bands,(/1,k%nibz/),title='-Fock Collisions')
 !
 isc%qs    =1
 isc%ngrho =N_g_vecs
 iscp%ngrho=N_g_vecs
 !
 ! Allocation
 !
 call elemental_collision_alloc(isc,NG=N_g_vecs,NG_GAMP=(/N_g_vecs,1/),TITLE="Fock")
 call elemental_collision_alloc(iscp,NG=N_g_vecs,TITLE="Fock")
 !
 call scatterGamp(isc,'x')
 !
 ! Parallel Indexes & LIVE-TIMING STEPS
 !
 do i_n = NL_bands(1),NL_bands(2)
   do i_m = i_n, NL_bands(2)
     !
     i_sp    = 1
     !
     Fock_collisions=cZERO
     !
     do iqbz=1,q%nbz
       !
       iqibz    =q%sstar(iqbz,1)
       iqs      =q%sstar(iqbz,2)
       !
       isc%qs   =(/0,iqibz,iqs/)
       iscp%qs  =(/0,iqibz,iqs/)
       !
       if (isc%qs(2)/=isc%iqref) call scatterGamp(isc,'x')
       !
       ! Scattering geometry (I)
       !=========================
       !
       !       (n,k,sp).     
       !               | m
       !               |   (q,r,sp)
       !               |
       !               | m'
       !    (n',kmq,sp).     
       !
       ! n   =QP_table(i_qp,1)
       ! m   =QP_table(i_qp,2)
       ! k   =QP_table(i_qp,3)
       !
       i_k_s   = 1
       !
       i_kmq_bz =qindx_S(i_k,iqbz,1)
       i_kmq_G0 =qindx_S(i_k,iqbz,2)
       !
       i_kmq   = k%sstar(i_kmq_bz,1)
       i_kmq_:qs = k%sstar(i_kmq_bz,2)
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
         call scatterBamp(iscp)
         !
         do i_np=n_bands(1),n_bands(2)
           !
           ! COLLISIONS(n,n',qp_n_states) = \sum_{G,q} ( rho(m,n,G,q,k) x rho*(m',n',G,q,k) /|q+G|^2 )
           !
           isc%os(1)=i_np
           !
           if (any(isc%is/=iscp%is).or.any(isc%os/=iscp%os)) then
             call scatterBamp(isc)
           else
             isc%rhotw=iscp%rhotw
           endif
           !
           pre_factor=cZERO
           !
           ! Remember isc%gamp=iscp%gamp
           !
!$omp parallel do default(shared), private(ig1), reduction(+:pre_factor)
           do ig1=1,N_g_vecs
             pre_factor=pre_factor+isc%rhotw(ig1)*isc%gamp(ig1,1)*conjg(iscp%rhotw(ig1))
           enddo
!$omp end parallel do
           !
           Fock_collisions(i_n,i_m,i_np,i_mp,i_kmq)=Fock_collisions(i_n,i_m,i_np,i_mp,i_kmq)-4._SP/spin_occ*pi*pre_factor
           !
         enddo ! loop on the band i_np
         !
       enddo ! loop on the band i_mp
       !
     enddo ! loop on iqbz
     !
     !
   enddo ! loop on i_m
 enddo ! loop on i_n
 !
 ! CLEAN
 !
 call collision_reset(isc)
 call collision_reset(iscp)
 !
end subroutine
