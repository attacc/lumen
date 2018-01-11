!
!        Copyright (C) 2000-2016 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM
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
subroutine COLLISIONS_free()
 !
 use pars,          ONLY:IP,schlen
 use electrons,     ONLY:levels,spin
 use R_lattice,     ONLY:bz_samp
 use wave_func,     ONLY:WF
 use collision,     ONLY:COLLISIONS_IO_status,COLLISIONS_state,Fock_collisions
 use SC,            ONLY:SC_potential
 use memory_m,      ONLY:mem_est
 !
 implicit none
 !
 if(allocated(COLLISIONS_IO_status)) then
   deallocate(COLLISIONS_IO_status)
   call mem_est("COLLISIONS_IO_status")
   deallocate(COLLISIONS_state)
   call mem_est("COLLISIONS_state")
 endif
 !
 if (trim(SC_potential)=="FOCK".or.trim(SC_potential)=="COHSEX") then
   !
   call COLLISIONS_alloc_and_free("Fock","all coll","clean",1)
   !
 endif
 !
 if (trim(SC_potential)=="COHSEX") then
   !
   call COLLISIONS_alloc_and_free("COHSEX","all coll","clean",1)
   !
 endif
 !
 ! CLEAN
 !=======
 !
 call WF_free(WF)
 !
end subroutine
