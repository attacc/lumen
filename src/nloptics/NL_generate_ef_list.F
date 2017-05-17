!
! Copyright (C) 2007 C. Attaccalite
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
subroutine NL_generate_ef_list(Efield)
 !
 ! Generate the list of Electric Fields used
 ! in the different runs
 !
 use pars ,        ONLY:SP
 use nl_optics,    ONLY:n_frequencies,n_intensities,n_runs,EF_parms,n_runs,nl_er,NL_estep
 use com,          ONLY:error
 use fields,       ONLY:ext_field,Efield_strength
 !
 implicit none
 !
 type(ext_field), intent(out) :: Efield
 !
 ! Work Space
 !
 integer  :: i_f,i_i,i_r
 real(SP) :: frequency
 !
 n_runs=n_frequencies*n_intensities
 !
 allocate(EF_parms(n_runs))
 !
 i_r=1
 !
 do i_f=1,n_frequencies
   !
   frequency=NL_er(1)+real(i_f-1._SP,SP)*NL_estep
   !
   do i_i=1,n_intensities
     !
     if(i_i==1) then
       EF_parms(i_r)%intensity =Efield%intensity
     elseif(i_i==2) then
       EF_parms(i_r)%intensity =Efield%intensity/4._SP
     elseif(i_i==3) then
       EF_parms(i_r)%intensity =Efield%intensity/16._SP
     else
       call error("Wrong intensity range! ")
     endif
     !
     EF_parms(i_r)%frequency = frequency
     EF_parms(i_r)%amplitude =Efield_strength(EF_parms(i_r)%intensity,'AU')
     EF_parms(i_r)%indexes   =(/i_f, i_i/)
     !
     i_r=i_r+1
     !
   enddo
   !
 enddo
 !
end subroutine NL_generate_ef_list

