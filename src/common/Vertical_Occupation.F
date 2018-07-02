!
!        Copyright (C) 2000-2017 Claudio Attaccalite
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
subroutine Vertical_Occupation(E,K,mode,ik_exc)
 !
 use pars,          ONLY:SP,DP
 use units,         ONLY:HA2EV
 use electrons,     ONLY:levels,nel_cond,nel,E_duplicate,E_reset
 use R_lattice,     ONLY:bz_samp
 use IO_m,          ONLY:IO_and_Messaging_switch
 use com,           ONLY:error,msg
 !
 implicit none
 !
 type(levels) ::E
 type(bz_samp)::K
 integer      ::mode
 integer, intent(in) :: ik_exc
 !
 ! Work Space
 !
 real(SP)     :: E_range=0.026/HA2EV !eV = 300 Kelvin
 integer :: n_bnd_hole,n_bnd_elect
 integer :: ib
 !
 do ib=1,E%nbf
   write(*,*) ' Energy ',E%E(ib,ik_exc,1)*HA2EV
 enddo
 do ib=E%nbf+1,E%nbf+10
   write(*,*) ' Energy ',E%E(ib,ik_exc,1)*HA2EV
 enddo
 stop
 !
end subroutine
