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
subroutine NL_long_range_SEX(dG,H,k,q,i_k) 
 !
 ! This subroutine calculate the long-range part of the Screened exchange
 !
 use electrons,      ONLY:n_sp_pol
 use pars,           ONLY:SP,cI,cZERO
 use FFT_m,          ONLY:fft_size
 use nl_optics,      ONLY:NL_bands,OSCL,W,COLL,COLL2
 use QP_m,           ONLY:QP_nk
 use R_lattice,      ONLY:bz_samp,qindx_S
 ! 
 implicit none
 !
 complex(SP), intent(inout) :: H(NL_bands(2),NL_bands(2),n_sp_pol)
 integer,     intent(in)    :: i_k
 complex(SP), intent(in)    :: dG(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),QP_nk)
 type(bz_samp), intent(in)  :: k,q
 !
 ! Work Space
 !
 integer :: i_m,i_mp  ! self-energy external indexes
 integer :: i_n,i_np,iqbz,i_kmq ! self-energy internal indexes
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
        do i_mp=NL_bands(1),NL_bands(2)
          do i_np=NL_bands(1),NL_bands(2)
            !
            H(i_n,i_m,1)=H(i_n,i_m,1)-cI*OSCL(i_n,i_np,i_k,iqbz)*conjg(OSCL(i_m,i_mp,i_k,iqbz))*W(iqbz) &
&                      *dG(i_np,i_mp,i_kmq) 
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
 enddo
 !
end subroutine NL_long_range_SEX
