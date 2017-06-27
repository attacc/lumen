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
subroutine OCCUPATIONS_Quasi_Fermi(E,K,mode)
 !
 ! Mode = 0 
 !  Update the Fermi energy
 !  Output : E%E_fermi
 !
 ! Mode = 1 -> Mode 0 +
 !  Define E%f and reports variations in Ef/nbf/nbm
 !  Output : E%nbf,E%nbm
 !  E%E are shifted 
 !
 ! Mode > 1 -> Mode 1 +
 !  Full report of system electronic character
 !
 use pars,          ONLY:SP,DP
 use electrons,     ONLY:levels,nel_cond
 use R_lattice,     ONLY:bz_samp
 implicit none
 !
 type(levels) ::E
 type(bz_samp)::K
 integer      ::mode
 !
 if(nel_cond==0._SP) then
   call  OCCUPATIONS_Fermi(E,K,mode)
   return
 endif
 !
end subroutine
