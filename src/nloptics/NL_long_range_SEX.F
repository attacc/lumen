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
subroutine NL_long_range_SEX(V_bands,H,k,q,E,i_k) 
 !
 ! This subroutine calculate the long-range part of the Screened exchange
 !
 use electrons,      ONLY:n_sp_pol,levels,spin_occ
 use pars,           ONLY:SP,cI,cZERO,cONE
 use FFT_m,          ONLY:fft_size
 use nl_optics,      ONLY:NL_bands,OSCL,W,ng_oscll,NL_nbands
 use QP_m,           ONLY:QP_nk
 use R_lattice,      ONLY:bz_samp,qindx_S
 use wrapper_omp,    ONLY:M_by_M_omp,Vstar_dot_V_omp
 ! 
 implicit none
 !
 complex(SP), intent(inout) :: H(NL_bands(2),NL_bands(2),n_sp_pol)
 integer,     intent(in)    :: i_k
 type(bz_samp), intent(in)  :: k,q
 type(levels),  intent(in)  :: E
 complex(SP),  intent(in)   :: V_bands(NL_bands(2),E%nbf,QP_nk,n_sp_pol)
 !
 ! Work Space
 !
 integer :: i_m,i_mp  ! self-energy external indexes
 integer :: i_n,i_np,iqbz,i_kmq ! self-energy internal indexes
 integer :: ig1,ig2
 complex(SP) :: COLL_element,dGelement
 complex(SP) OSCL_val(NL_bands(1):NL_bands(2),E%nbf,q%nbz)
 !
 !
 ! Build Oscillators in term of valence bands
 ! 
 do i_n=NL_bands(1),NL_bands(2)
   do iqbz=1,q%nbz
     i_kmq   = k%sstar(qindx_S(i_k,iqbz,1),1)
     call M_by_M_omp('C','N',NL_nbands,E%nbf,NL_nbands,cONE,OSCL(1,:,:,i_k,iqbz),&
&                    NL_nbands,V_bands(:,:,i_kmq,1),NL_bands(2),cZERO,OSCL_val(:,:,iqbz),NL_nbands)

   enddo
 enddo
 !
 ! Long-range SEX
 !
 !
 do i_n=NL_bands(1),NL_bands(2)
   do i_m=i_n,NL_bands(2)
     !
     ! Evaluate Sigma_SEX(i_m,i_n,i_k)
     !
     do iqbz=1,q%nbz
        i_kmq   = k%sstar(qindx_S(i_k,iqbz,1),1)
        !
!        call M_by_M_omp('C','N',NL_bands(2),en%nbf,NL_bands(2),cONE,OSCL(1,:,:,i_k,iqbz),&
!&                          NL_bands(2),V_bands(:,:,i_kmq),NL_bands(2),cZERO,OCSL_val,NL_bands(2))
        
        do i_mp=NL_bands(1),NL_bands(2)
          do i_np=NL_bands(1),NL_bands(2)
            !
            ! Build Collision on the fly
            !
            COLL_element=cZERO
            do ig1=1,ng_oscll
              do ig2=1,ng_oscll
                COLL_element=COLL_element+OSCL(ig1,i_n,i_np,i_k,iqbz)*conjg(OSCL(ig2,i_m,i_mp,i_k,iqbz))*W(iqbz,ig2,ig1) 
              enddo
            enddo
            !
            dGelement=cI*spin_occ*Vstar_dot_V_omp(E%nbf,V_bands(i_mp,:,i_kmq,1),V_bands(i_np,:,i_kmq,1))
            if(i_mp==i_np.and.i_np<=E%nbf) dGelement=dGelement-cI*spin_occ
            H(i_n,i_m,1)=H(i_n,i_m,1)-cI*COLL_element*dGelement 
            ! 
          enddo
        enddo
        !
     enddo
     !
   enddo
 enddo
 !
 ! Symmetrization
 !
 do i_n=NL_bands(1),NL_bands(2)
   do i_m=i_n+1,NL_bands(2)
     H(i_m,i_n,1)=conjg(H(i_n,i_m,1))
   enddo
   H(i_n,i_n,1)=real(H(i_n,i_n,1))
 enddo
 !
end subroutine NL_long_range_SEX
