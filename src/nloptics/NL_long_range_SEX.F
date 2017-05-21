!
!        Copyright (C) 2000-2015 the LUMEN team
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
subroutine NL_long_range_SEX(V_bands,Hk,k,q,E,i_k) 
 !
 ! This subroutine calculate the long-range part of the Screened exchange
 !
 use electrons,      ONLY:n_sp_pol,levels,spin_occ
 use pars,           ONLY:SP,cI,cZERO,cONE
 use FFT_m,          ONLY:fft_size
 use nl_optics,      ONLY:NL_bands,OSCLL,W,ng_oscll,NL_nbands,Sigma_SEX_EQ
 use QP_m,           ONLY:QP_nk
 use R_lattice,      ONLY:bz_samp,qindx_S
 use wrapper,        ONLY:M_by_M,Vstar_dot_V,M_by_V
 ! 
 implicit none
 !
 complex(SP), intent(inout) :: Hk(NL_bands(2),NL_bands(2),n_sp_pol)
 integer,     intent(in)    :: i_k
 type(bz_samp), intent(in)  :: k,q
 type(levels),  intent(in)  :: E
 complex(SP),  intent(in)   :: V_bands(NL_bands(2),E%nbf,QP_nk,n_sp_pol)
 !
 ! Work Space
 !
 integer :: i_m,i_mp,i_v  ! self-energy external indexes
 integer :: i_n,i_np,iqbz,i_kmq ! self-energy internal indexes
 integer :: ig1,ig2,i1
 complex(SP) :: OSCLL_val(NL_bands(1):NL_bands(2),E%nbf,ng_oscll,q%nbz)
 complex(SP) :: qGG_sum(E%nbf)
 !
 !
 ! Build Oscillators in term of valence bands
 ! 
!$omp parallel do default(shared), private(iqbz,i_kmq,ig1)
 do iqbz=1,q%nbz
   i_kmq   = k%sstar(qindx_S(i_k,iqbz,1),1)
   do ig1=1,ng_oscll
     call M_by_M('N','N',NL_nbands,E%nbf,NL_nbands,cONE,OSCLL(:,:,i_k,iqbz,ig1),&
&                 NL_nbands,V_bands(NL_bands(1):,:,i_kmq,1),NL_nbands,cZERO,OSCLL_val(:,:,ig1,iqbz),NL_nbands)
   enddo
 enddo
!$omp end parallel do
 !
!$omp parallel do default(shared), private(i1)
 do i1=1,NL_bands(2)
   Hk(:,i1,:)=Hk(:,i1,:)-Sigma_SEX_EQ(:,i1,:,i_k)
 enddo
!$omp end parallel do
 !
 ! Long-range SEX
 !
 ! Evaluate Sigma_SEX(i_m,i_n,i_k)
 !
 do iqbz=1,q%nbz
   i_kmq   = k%sstar(qindx_S(i_k,iqbz,1),1)
   !
   do ig1=1,ng_oscll
     do ig2=1,ng_oscll
       call M_by_M('N','C',NL_nbands,NL_nbands,E%nbf,spin_occ*W(iqbz,ig1,ig2),OSCLL_val(:,:,ig1,iqbz),&
&                      NL_nbands,OSCLL_val(:,:,ig2,iqbz),NL_nbands,cONE,Hk(NL_bands(1):,NL_bands(1):,1),NL_nbands)
     enddo
   enddo
   !
 enddo
 ! 
! do i_n=NL_bands(1),NL_bands(2)
!   do i_m=i_n,NL_bands(2)
!     !
!     do i_v=1,E%nbf
!       do ig1=1,ng_oscll
!         do ig2=1,ng_oscll
!          qGG_sum(i_v)=sum(W(:,ig1,ig2)*OSCLL_val(i_n,i_v,ig1,:)*conjg(OSCLL_val(i_m,i_v,ig2,:)))
!         enddo
!       enddo
!     enddo
!     !
!     Hk(i_n,i_m,1)=Hk(i_n,i_m,1)+spin_occ*sum(qGG_sum)
!     !
!   enddo
! enddo
! !
! ! Symmetrization
! !
! do i_n=NL_bands(1),NL_bands(2)
!   do i_m=i_n+1,NL_bands(2)
!     Hk(i_m,i_n,1)=conjg(Hk(i_n,i_m,1))
!   enddo
!   Hk(i_n,i_n,1)=real(Hk(i_n,i_n,1))
! enddo
 !
end subroutine NL_long_range_SEX
