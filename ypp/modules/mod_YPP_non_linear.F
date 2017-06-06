!
!        Copyright (C) 2000-2016 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AC
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
module YPP_non_linear
 !
 use pars,     ONLY:schlen,SP,DP
 use units,    ONLY:FS2AUT,HA2EV
 !
 implicit none
 !
 ! (NL) optics
 !
 logical  :: Use_FFT
 logical  :: l_force_min_damp
 integer  :: N_probe_frequencies
 integer  :: X_order
 real(SP) :: W_minimal_step ! minimal frequency sted defined from the simulation length
 integer  :: N_richardson        =1  ! order of Richardson extrapolation
 !
 type NL_resp
   real(SP) :: delta_T           =0._SP                    ! delta T
   real(SP) :: T_range(2)        =-FS2AUT                  ! time windows for the transform
   integer  :: i_T_range(2)      =-1                       ! integer time windows 
   real(SP) :: W_range(2)        =(/0._SP,10._SP/HA2EV/)   ! energy range 
   integer  :: T_n_steps         =0                        ! number of time steps for the response function
   integer  :: W_n_steps         =200                      ! number of steps in energy mesh
   real(SP) :: damp_factor       =0.1_SP/HA2EV   ! factor used in damping
   real(SP) :: Field_damp_factor =0._SP   ! factor used deduced from the Field 
   character(schlen) :: damping  ="NONE"  ! Damp type (none, lorentian or gaussian)
 end type NL_resp
 !
 type(NL_resp), save :: NL_response
 !
 complex(SP), allocatable :: Div_Efield(:)
 !
end module YPP_non_linear
