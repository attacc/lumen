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
subroutine NL_lsex_potential(q) 
 !
 ! Construct the screeend potential W(G=G'=0,q)
 !
 use nl_optics,           ONLY:NL_bands,W
 use LIVE_t,              ONLY:live_timing
 use pars,                ONLY:SP,cZERO
 use frequency,           ONLY:w_samp,W_reset
 use R_lattice,           ONLY:bz_samp,qindx_S
 use electrons,           ONLY:levels
 use collision,           ONLY:elemental_collision,elemental_collision_free, &
&                              elemental_collision_alloc
 use X_m,                 ONLY:X_mat,X_alloc,X_t
 use IO_m,                ONLY:io_control,OP_RD_CL,REP,VERIFY,NONE,RD_CL,OP_RD,RD_CL_IF_END 
 use com,                 ONLY:error
 !
 implicit none
 !
 type(bz_samp), intent(in) :: q
 !
 ! Work Space
 !
 integer                    :: iqbz,iqibz,iqs,io_err,ID
 integer, external          :: io_X
 type(elemental_collision)  :: isc 
 type(w_samp)               :: Xw  
 type(X_t)                  :: X
 !
 call W_reset(Xw)
 Xw%n_freqs=1
 call FREQUENCIES_setup(Xw)
 !
 X%whoami=2 ! load em1s
 !
 isc%iqref =0
 isc%ngrho =1  ! only G=1 is required for local-SEX 
 !
 call elemental_collision_free(isc)
 call elemental_collision_alloc(isc,NG_GAMP=(/1,1/),TITLE="Wlr") 
 !
 call X_alloc('X',(/1,1,1/))
 ! 
 call live_timing('[NL] W potential :',q%nbz)
 !
 do iqbz=1,q%nbz
   ! 
   isc%qs(2:)=(/q%sstar(iqbz,1),q%sstar(iqbz,2)/)
   iqibz=isc%qs(2)
   iqs  =isc%qs(3)
   !
   call scatter_Gamp(isc,'c')
   !
   if (iqbz ==1) call io_control(ACTION=OP_RD,COM=NONE,       SEC=(/1,2,2*iqibz+1/),ID=ID)
   if (q%nbz==1) call io_control(ACTION=OP_RD_CL,COM=NONE,    SEC=(/1,2,3/),ID=ID)
   if (iqbz>  1) call io_control(ACTION=RD_CL_IF_END,COM=NONE,SEC=(/2*iqibz,2*iqibz+1/),ID=ID)
   io_err=io_X(X,Xw,ID)
   !
   if(io_err/=0) call error("Error loading em1s")
   !
   W(iqbz)=X_mat(1,1,1)*isc%gamp(1,1)
   !
   call live_timing(steps=1)
   !
 enddo
 !
 call live_timing()
 !
 call elemental_collision_free(isc)
 call X_alloc('X')
 !
end subroutine NL_lsex_potential
