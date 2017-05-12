!
!        Copyright (C) 2000-2016 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM CA
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
subroutine NL_Rdamp_it(damp_type,damp_factor,ft,lda)
 !
 use YPP_non_linear, ONLY:NL_response
 use pars,           ONLY:SP,DP
 use com,            ONLY:warning
 use nl_optics,      ONLY:NL_step

 implicit none
 !
 integer,      intent(in)      :: lda
 real(SP),     intent(in)      :: damp_factor
 real(DP),     intent(inout)   :: ft(NL_response%T_n_steps,lda)
 character(*), intent(in)      :: damp_type
 !                          
 ! Work Space 
 !
 integer  :: i_t,i_l
 real(SP) :: damp_function(NL_response%T_n_steps)
 !
 do i_t=1,NL_response%T_n_steps
   select case ( trim(damp_type) )
     case("LORENTZIAN")  
       damp_function(i_t) = exp(-abs(NL_step*real(i_t,SP)*damp_factor))
     case("GAUSSIAN") 
       damp_function(i_t) = exp(-(NL_step*real(i_t,SP))**2*damp_factor**2)
     case("NONE") 
       damp_function(i_t) = 1._SP
     case default
       call warning(" Unknown damping! ")
       return
   end select
 enddo
 !
 forall(i_l=1:lda)  ft(:,i_l)=ft(:,i_l)*damp_function(:)
 !
end subroutine NL_Rdamp_it

