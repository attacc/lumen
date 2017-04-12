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
subroutine NL_lsex_potential(q,X,Xw) 
 !
 ! Construct the screeend potential W(G=G'=0,q)
 !
 use nl_optics,           ONLY:NL_bands,W,Correlation,LSEX
 use LIVE_t,              ONLY:live_timing
 use pars,                ONLY:SP,cZERO,pi
 use frequency,           ONLY:w_samp,W_reset
 use R_lattice,           ONLY:bz_samp,qindx_S
 use D_lattice,           ONLY:DL_vol
 use electrons,           ONLY:levels,spin_occ
 use collision,           ONLY:elemental_collision,elemental_collision_free, &
&                              elemental_collision_alloc
 use X_m,                 ONLY:X_mat,X_alloc,X_t
 use IO_m,                ONLY:io_control,OP_RD_CL,REP,VERIFY,NONE,RD_CL,OP_RD,RD_CL_IF_END 
 use com,                 ONLY:error
 use stderr,              ONLY:intc
 !
 implicit none
 !
 type(bz_samp), intent(in) :: q
 type(w_samp)              :: Xw
 type(X_t)                 :: X
 !
 ! Work Space
 !
 integer                    :: iqbz,iqibz,iqs,io_err,ID
 integer, external          :: io_X
 type(elemental_collision)  :: isc 
 !
 isc%iqref =0
 isc%ngrho =X%ng  ! only G=1 is required for local-SEX 
 !
 call elemental_collision_free(isc)
 call elemental_collision_alloc(isc,NG_GAMP=(/X%ng,1/),TITLE="Wlr") 
 !
 call X_alloc('X',(/X%ng,X%ng,1/))
 ! 
 call live_timing('[NL] W potential :',q%nbz)
 !
 if(Correlation==LSEX) then
   call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1,2/),MODE=VERIFY,ID=ID)
   io_err=io_X(X,Xw,ID)
   if (io_err<0) call error('Incomplete and/or broken tatic diel. fun. database')
 endif
 !
 do iqbz=1,q%nbz
   ! 
   isc%qs(2:)=(/q%sstar(iqbz,1),q%sstar(iqbz,2)/)
   iqibz=isc%qs(2)
   iqs  =isc%qs(3)
   !
   isc%iqref=-1
   call scatter_Gamp(isc,'x')
   !
   ! The bare exchange v
   !
   W(iqbz)=-4._SP/spin_occ*pi*isc%gamp(1,1)
   !
   if(Correlation==LSEX) then
     !      
     isc%iqref=-1
     call scatter_Gamp(isc,'c')
     !
     if (iqbz ==1) then
       call io_control(ACTION=OP_RD,COM=NONE,       SEC=(/1,2,2*iqibz+1/),ID=ID)
     elseif (q%nbz==1) then
       call io_control(ACTION=OP_RD_CL,COM=NONE,    SEC=(/1,2,3/),ID=ID)
     elseif (iqbz>  1) then 
       call io_control(ACTION=RD_CL_IF_END,COM=NONE,SEC=(/2*iqibz,2*iqibz+1/),ID=ID)
     endif
     io_err=io_X(X,Xw,ID)
     !
     if(io_err<0) call error(" Error loading em1s iq="//intc(iqbz))
     !
     !  Screened part of the exchange  \tilde W = W - v
     !
     W(iqbz)=W(iqbz)-4._SP/spin_occ*pi*X_mat(1,1,1)*isc%gamp(1,1)
     !
   endif
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
