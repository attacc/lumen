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
!
subroutine NL_1D_Fourier_Transform(mode,title,freqs,fw,nw,ft,nt,I_sign)
 !
 ! Here I calculate:
 !
 ! mode="T2W"
 !
 ! \tilde{FT}[ft] = \int_{T0}{\infty} ft(t) exp^(i I_sign w t) =  FT[ft]   
 !
 ! This definition follows from the fact that ft(t) is defined in input only for t>0.
 !
 ! mode="W2T" 
 !
 ! ft(t) = \int fw(w) exp^(-i w t) 
 !
 use YPP_non_linear, ONLY:NL_response
 use nl_optics,      ONLY:NL_step
 use pars,           ONLY:SP,DP,cI
 use LIVE_t,         ONLY:live_timing
 implicit none
 !
 integer,     intent(in)   :: nw,nt,I_sign
 complex(SP), intent(in)   :: freqs(nw)
 character(*),intent(in)   :: mode
 character(*),intent(in)   :: title
 complex(DP), intent(inout):: ft(nt)
 complex(DP), intent(inout):: fw(nw)
 !
 ! Work Space 
 !
 integer    :: i_w,i_t
 logical    :: l_live_timing
 real(SP)   :: NL_time
 !
 l_live_timing=len_trim(title)>0
 !
 ! Time to Frequency
 !===================
 if (mode=="T2W") then
   !
   fw=(0._DP,0._DP)
   !
   if (l_live_timing) call live_timing('FT'//title//':',nw)
   do i_w = 1, nw
     do i_t = 1, nt 
        NL_time=(i_t-1)*NL_step
        fw(i_w) = fw(i_w)+ft(i_t)*exp(cI*I_sign*freqs(i_w)*NL_time)*NL_step
     enddo
     if (l_live_timing) call live_timing(steps=1)
   enddo
   if (l_live_timing) call live_timing()
   !
 endif
 !
 ! Frequency to Time
 !===================
 if (mode=="W2T") then
   !
   ft=(0._DP,0._DP)
   !
   if (l_live_timing) call live_timing('FT'//title//':',nt)
   do i_t = 1, nt
     do i_w = 1, nw
     !       ft(i_t) = ft(i_t)+fw(i_w)*(freqs(2)-freqs(1))*exp(-cI*freqs(i_w)*RT_time(i_t))
     enddo
     if (l_live_timing) call live_timing(steps=1)
   enddo
   if (l_live_timing) call live_timing()
   !
 endif
 !
end subroutine NL_1D_Fourier_Transform
