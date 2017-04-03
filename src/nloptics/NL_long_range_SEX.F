!
!        Copyright (C) 2000-2015 the LUMEN team
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
subroutine NL_long_range_SEX(dG,H,ik) 
 !
 ! This subroutine calculate the long-range part of the Screened exchange
 !
 use electrons,      ONLY:n_sp_pol
 use pars,           ONLY:SP
 use FFT_m,          ONLY:fft_size
 use nl_optics,      ONLY:NL_bands
 use QP_m,        ONLY:QP_nk
 !
 complex(SP), intent(inout) :: H(NL_bands(2),NL_bands(2),n_sp_pol)
 integer,     intent(in)    :: ik
 complex(SP), intent(in)    :: dG(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),QP_nk)
 !
end subroutine NL_long_range_SEX
