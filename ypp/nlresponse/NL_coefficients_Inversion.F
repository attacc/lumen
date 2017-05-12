!
!        Copyright (C) 2000-2016 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM MG
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
subroutine NL_coefficients_Inversion(NW,NX,P,X,W,T_period)
 !
 ! See reference
 ! PRB 88 235113 (2013) 
 !
 use pars,            ONLY:SP,cI,cONE,cZERO
 use YPP_non_linear,  ONLY:X_order,NL_response
 use nl_optics,       ONLY:NL_step,NL_steps
 use linear_algebra,  ONLY:LINEAR_ALGEBRA_driver,INV
 implicit none
 !
 integer,     intent(in)   :: NW,NX
 real(SP),    intent(in)   :: P(NL_steps),W(NW),T_period
 complex(SP), intent(out)  :: X(NX)
 !
 ! Work space
 !
 integer                 :: i_t,i_n,i_deltaT,nP_components,M_size
 real(SP)   ,allocatable :: P_i(:),T_i(:)
 complex(SP),allocatable :: M(:,:)
 !
 M_size        =2*X_order+1
 nP_components =X_order+1 
 !
 i_deltaT=nint(T_period/NL_step)/M_size
 !
 allocate(M(M_size,M_size),P_i(M_size),T_i(M_size))
 !
 M=cZERO
 !
 do i_t=1,M_size
   T_i(i_t)= ((NL_response%i_T_range(1)-1._SP)+i_deltaT*(i_t-1._SP))*NL_step
   P_i(i_t)=P(NL_response%i_T_range(1) + (i_t-1)*i_deltaT)
 enddo
 !
 !
 do i_t=1,M_size
   M(i_t,1) = cONE
   do i_n=2,nP_components
     M(i_t,i_n        ) = exp( cI*W(i_n)*T_i(i_t))
     M(i_t,i_n+X_order) = exp(-cI*W(i_n)*T_i(i_t))
   end do
 end do
 !
 call LINEAR_ALGEBRA_driver(INV,M)
 !
 X = cZERO
 do i_n =1, nP_components
   do i_t = 1,M_size
     X(i_n)=X(i_n)+M(i_n,i_t)*P_i(i_t) 
   end do
 end do
 !
 deallocate(M,P_i,T_i)
 ! 
end subroutine NL_coefficients_Inversion

