!
!        Copyright (C) 2000-2017 the YAMBO team
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
subroutine NL_oscillators_load(q) 
 !
 ! Evaluate oscillators for the Local-SEX
 !
 use LIVE_t,              ONLY:live_timing
 use R_lattice,           ONLY:bz_samp
 use QP_m,                ONLY:QP_nk
 use parallel_m,          ONLY:PAR_IND_Xk_ibz,PAR_IND_Xk_ibz_ID
 use IO_m,                ONLY:io_control,OP_RD_CL,VERIFY,REP
 use com,                 ONLY:msg
 use parallel_m,          ONLY:PP_redux_wait,PAR_IND_Xk_ibz,PAR_COM_Xk_ibz_INDEX,master_cpu
 use IO_m,                ONLY:IO_and_Messaging_switch
 !
 implicit none
 !
 type(bz_samp), intent(in) :: q
 !
 ! Work Space
 !
 integer :: i_k
 integer :: ID,io_err
 integer, external :: io_OSCLL
 !
 call PP_redux_wait() 
 !
 i_sp_pol=1
 !
 call IO_and_Messaging_switch("SAVE")
 call IO_and_Messaging_switch("+io_in",CONDITION=.TRUE.)
 !   
 call msg('s','[OSC] Checking oscillators header ')
 call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1/),MODE=VERIFY,ID=ID)
 io_err=io_OSCLL(q,ID)
 ! 
 ! Load Oscillators
 !
 call live_timing('[OSC] Load Oscillators :',PAR_IND_Xk_ibz%n_of_elements(PAR_IND_Xk_ibz_ID+1))
 !
 do i_k=1,QP_nk
   if (.not.PAR_IND_Xk_ibz%element_1D(i_k)) cycle
   call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1+i_k+(i_sp_pol-1)*k%nibz/),ID=ID)
   io_err=io_OSCLL(q,ID)
   call live_timing(steps=1)
 enddo
 call live_timing()
 !
 call IO_and_Messaging_switch("RESTORE")
 !
 call msg('rs','[OSC] Oscillators loaded')
 !
end subroutine NL_oscillators_load
