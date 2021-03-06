! 
!        Copyright (C) 2000-2016 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): DS
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
subroutine electrons_spin_and_magn_interpolate(Xk,CIRCUIT_k,BANDS_range,spin_output,magn_output)
 !
 use pars,           ONLY:SP,schlen
 use units,          ONLY:HA2EV
 use electrons,      ONLY:levels,n_spinor
 use R_lattice,      ONLY:bz_samp,bz_samp_reset
 use com,            ONLY:msg,of_open_close,error,warning
 use vec_operate,    ONLY:c2a,v_norm,iku_v_norm
 use parser_m,       ONLY:parser
 use stderr,         ONLY:intc
 use interpolate,    ONLY:eval_interpolation_coeff,bz_interpolation
 use QP_CTL_m,       ONLY:QP_apply
 !
 implicit none
 !
 type(bz_samp), intent(in)  :: Xk
 integer      , intent(in)  :: BANDS_range(2)
 type(bz_samp), intent(in)  :: CIRCUIT_k
 real(SP)     , intent(out) :: spin_output(BANDS_range(1):BANDS_range(2),CIRCUIT_k%nibz,n_spinor),&
&                              magn_output(BANDS_range(1):BANDS_range(2),CIRCUIT_k%nibz,3)
 !
 ! Work space
 !
 integer  :: i1,ib,ik,i_sp,i_dir,n_states, &
&            tmp_table((BANDS_range(2)-BANDS_range(1)+1)*Xk%nibz,4)
 real(SP) :: spin_fac(BANDS_range(1):BANDS_range(2),Xk%nibz,n_spinor),        &
&            magn_fac(BANDS_range(1):BANDS_range(2),Xk%nibz,3),               &
&            tmp_factor_in((BANDS_range(2)-BANDS_range(1)+1)*Xk%nibz,5),      &
&            tmp_factor_out(BANDS_range(1):BANDS_range(2),CIRCUIT_k%nibz)
 !
 call electrons_spin_factors(Xk,(/1,Xk%nibz/),BANDS_range,spin_fac)
 !
 call electrons_magn_factors(Xk,(/1,Xk%nibz/),BANDS_range,magn_fac)
 !
 n_states=(BANDS_range(2)-BANDS_range(1)+1)*Xk%nibz
 tmp_table=0
 !
 i1=0
 do ik=1,Xk%nibz
   do ib=BANDS_range(1),BANDS_range(2)
     i1=i1+1
     tmp_table(i1,1:2)=ib
     tmp_table(i1,3)  =ik
     tmp_factor_in(i1,1:2)=spin_fac(ib,ik,:)
     tmp_factor_in(i1,3:5)=magn_fac(ib,ik,:)
   enddo
 enddo
 ! 
 do i_sp=1,n_spinor
   tmp_factor_out=0._SP
   call Nearest_kpt_interpolation(Xk%nibz,CIRCUIT_k%nibz,BANDS_range,n_states,1,1,tmp_table,Xk%pt,CIRCUIT_k%pt, &
&                                   Xk%weights,CIRCUIT_k%weights,tmp_factor_in(:,i_sp),tmp_factor_out,'SPN',.false.)
   do ik=1,CIRCUIT_k%nibz
     spin_output(:,ik,i_sp)=tmp_factor_out(:,ik)
   enddo
 enddo
 !
 do i_dir=1,3
   tmp_factor_out=0._SP
   call Nearest_kpt_interpolation(Xk%nibz,CIRCUIT_k%nibz,BANDS_range,n_states,1,1,tmp_table,Xk%pt,CIRCUIT_k%pt, &
&                                   Xk%weights,CIRCUIT_k%weights,tmp_factor_in(:,2+i_dir),tmp_factor_out,'MAG',.false.)
   do ik=1,CIRCUIT_k%nibz
     magn_output(:,ik,i_dir)=tmp_factor_out(:,ik)
   enddo
 enddo
 !
end subroutine electrons_spin_and_magn_interpolate
