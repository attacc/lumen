!
!        Copyright (C) 2000-2016 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM CA
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
subroutine V_real_space_to_H(ik,V,H,WF,V_kind)
 !
 use electrons,      ONLY:n_sp_pol,n_spin,n_spinor,spin
 use pars,           ONLY:SP
 use FFT_m,          ONLY:fft_size
 use wave_func,      ONLY:WAVEs
 use wrapper_omp,    ONLY:Vstar_dot_V_omp
 use matrix_operate, ONLY:UP_matrix_index
 use nl_optics,      ONLY:NL_bands
 implicit none
 !
 integer     :: ik,ir,i_spin
 real(SP)    :: V(fft_size,n_spin)
 complex(SP) :: V_mat(fft_size,n_spin,n_spin)
 complex(SP) :: H(NL_bands(2),NL_bands(2),n_sp_pol)
 character(3):: V_kind
 type(WAVEs) :: WF
 !
 ! Work Space
 !
 integer    :: i_qp,ib,ibp,ifft,ifftp
 complex(SP):: c_ws(fft_size,n_spin)
 !
 if(V_kind=="xc ") call Build_V_xc_mat(V,V_mat)
 !
 ! Fill ...
 !
 do i_spin=1,n_sp_pol
   do ib=NL_bands(1),NL_bands(2)
     do ibp=ib,NL_bands(2)
     !
     !
     ! <ifft|V|ifftp>
     !
     if (n_spinor==1) then
       !
       ifft   =WF%index(ib,ik,i_spin)
       ifftp  =WF%index(ibp,ik,i_spin)
!$omp parallel do default(shared), private(ir)
       do ir=1,fft_size
         c_ws(ir,1)=V(ir,i_spin)*WF%c(ir,ifftp)
       enddo
!$omp end parallel do
       !
       H(ib,ibp,i_spin) = H(ib,ibp,i_spin)+Vstar_dot_V_omp(fft_size,WF%c(:,ifft),c_ws(:,1))
       !
     elseif (n_spinor==2) then
       !
       if (V_kind=='def') then
         !
         ifft   =WF%index(ib,ik,1)
         ifftp  =WF%index(ibp,ik,1)
!$omp parallel do default(shared), private(ir)
         do ir=1,fft_size
           c_ws(ir,1)=V(ir,1)*WF%c(ir,ifftp)
         enddo
!$omp end parallel do
         H(ib,ibp,1) = H(ib,ibp,1)+Vstar_dot_V_omp(fft_size,WF%c(:,ifft),c_ws(:,1))
         ifft   =WF%index(ib,ik,2)
         ifftp  =WF%index(ibp,ik,2)
!$omp parallel do default(shared), private(ir)
         do ir=1,fft_size
           c_ws(ir,2)=V(ir,2)*WF%c(ir,ifftp)
         enddo
!$omp end parallel do
         H(ib,ibp,1) = H(ib,ibp,1)+Vstar_dot_V_omp(fft_size,WF%c(:,ifft),c_ws(:,2))
         !
       elseif(V_kind=='xc ') then
         !
         ifft   =WF%index(ib,ik,1)
         !
         ifftp  =WF%index(ibp,ik,1)
!$omp parallel do default(shared), private(ir)
         do ir=1,fft_size
           c_ws(ir,1)=V_mat(ir,1,1)*WF%c(ir,ifftp)
         enddo
!$omp end parallel do
         ifftp  =WF%index(ibp,ik,2)
!$omp parallel do default(shared), private(ir)
         do ir=1,fft_size
           c_ws(ir,2)=V_mat(ir,1,2)*WF%c(ir,ifftp)
         enddo
!$omp end parallel do
         H(ib,ibp,1) = H(ib,ibp,1)+Vstar_dot_V_omp(fft_size,WF%c(:,ifft),c_ws(:,1))+&
&                                  Vstar_dot_V_omp(fft_size,WF%c(:,ifft),c_ws(:,2))
         !
         ifft   =WF%index(ib,ik,2)
         !
         ifftp  =WF%index(ibp,ik,1)
!$omp parallel do default(shared), private(ir)
         do ir=1,fft_size
           c_ws(ir,1)=V_mat(ir,2,1)*WF%c(ir,ifftp)
         enddo
!$omp end parallel do
         ifftp  =WF%index(ibp,ik,2)
!$omp parallel do default(shared), private(ir)
         do ir=1,fft_size
           c_ws(ir,2)=V_mat(ir,2,2)*WF%c(ir,ifftp)
         enddo
!$omp end parallel do
         H(ib,ibp,1) = H(ib,ibp,1)+Vstar_dot_V_omp(fft_size,WF%c(:,ifft),c_ws(:,1))+&
&                                  Vstar_dot_V_omp(fft_size,WF%c(:,ifft),c_ws(:,2))
         !
       endif
       !
     endif
     !
   enddo
  enddo
 enddo
 !
 ! ... and symmetrize
 !
!$omp parallel do default(shared), private(ib,ibp)
 do ib=NL_bands(1),NL_bands(2)
   H(ib,ib,:)=real(H(ib,ib,:),SP)
   do ibp=1,ib-1
     H(ib,ibp,:)=conjg(H(ibp,ib,:))
   enddo
 enddo
!$omp end parallel do
 !
end subroutine 
