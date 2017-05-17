!
!        Copyright (C) 2007 Claudio Attaccalite and Myrta Gruning
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
subroutine NL_richardson(Susc_vs_field,Susc_extrapolated)
 !
 use pars,           ONLY:SP
 use nl_optics,      ONLY:n_runs,n_frequencies,n_intensities
 use YPP_non_linear, ONLY:X_order
 use com,            ONLY:error
 !
 implicit none
 !
 complex(SP), intent(in)  :: Susc_vs_field(0:X_order,3,n_frequencies,n_intensities)
 complex(SP), intent(out) :: Susc_extrapolated(0:X_order,3,n_frequencies)
 !
 if(n_intensities==1) then
   !
   ! No richardson extrapolation
   !   
   Susc_extrapolated=Susc_vs_field(:,:,:,1)
   !
 elseif(n_intensities==2) then
   !
   ! For a system with inversion symmetry
   ! Remove linear component 
   !
   ! X^(2)(2w) = ..... work in progess 
   ! X^(3)(3w) = [ 4 * P(3w,E)/E^(3) - P(3w, E/2)/ (E/2)^3 ] /3 + O(5)
   !
   Susc_extrapolated=(4._SP*Susc_vs_field(:,:,:,1) - Susc_vs_field(:,:,:,2))/3._SP
   !
 elseif(n_intensities==3) then
   !
   ! For a system without inversion symmetry 
   ! Remove linear and quadratic components 
   !
   ! P(w,E)   = a E   + b E^2    + X3 E^3
   ! P(w,E/2) = a E/2 + b E^2/4  + X3 E^3/8
   ! P(w,E/4) = a E/4 + b E^2/16 + X3 E^3/16
   !
   ! X3 = 8/3 [ P(w,E) - 6 P(w, E/2) + 8 P(w, E/4) ]/E^3
   !
   Susc_extrapolated=Susc_vs_field(:,:,:,1) - 6._SP*Susc_vs_field(:,:,:,2) + 8._SP*Susc_vs_field(:,:,:,3)
   Susc_extrapolated=8._SP/3._SP*Susc_extrapolated
   !
 else
   call error("Richardson extrapolation order==4 not implemented!! ")
 endif
 !
end subroutine NL_richardson
