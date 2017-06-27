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
subroutine REPORT_all_Occupations(E,k,what)
 !
 use pars,           ONLY:schlen,SP
 use units,          ONLY:HA2EV
 use R_lattice,      ONLY:bz_samp,nXkibz
 use electrons,      ONLY:n_sp_pol,nel,levels
 use stderr,         ONLY:intc,real2ch
 use wave_func,      ONLY:wf_nc_k
 use com,            ONLY:msg
 use vec_operate,    ONLY:c2a
 !
 implicit none
 !
 type(bz_samp) :: k
 type(levels)  :: E
 character(*)  :: what
 !
 ! Work Space
 !
 integer           :: ik,is_pol,ib,i1
 character(schlen) :: dumb_ch1,dumb_ch2
 real(SP)      :: k_cc(3),k_al(3)
 !
 call section('Occupations Report')
 !
 do ik=1,E%nk
   !
   call c2a(v_in=k%pt(ik,:),v_out=k_al,mode='ki2a')
   call c2a(v_in=k%pt(ik,:),v_out=k_cc,mode='ki2c')
   ! 
   dumb_ch1='K ['//trim(intc(ik))//'] :'//trim(real2ch(k_cc(1)))//' '//&
&          trim(real2ch(k_cc(2)))//' '//trim(real2ch(k_cc(3)))//&
&          ' ( cc) * Comp.s '//trim(intc(wf_nc_k(ik)))
   dumb_ch2='       '//trim(real2ch(k%pt(ik,1)))//' '//&
&          trim(real2ch(k%pt(ik,2)))//' '//trim(real2ch(k%pt(ik,3)))//&
&          ' (rlu) '
   if(ik>  9) dumb_ch2=' '//trim(dumb_ch2)
   if(ik> 99) dumb_ch2=' '//trim(dumb_ch2)
   if(ik>999) dumb_ch2=' '//trim(dumb_ch2)
   !
   if (ik> 1    ) call msg('r','')
   call msg('r','*X* '//trim(dumb_ch1)//' * weight '//trim(real2ch(k%weights(ik))))
   call msg('r','    '//trim(dumb_ch2))
   !
   do is_pol=1,n_sp_pol
     !
     do i1=1,E%nb,8
       if (n_sp_pol==1) then
         call msg('r',what//' ',(/(E%E(ib,ik,1),ib=i1,min(i1+7,E%nb))/))
       else
         if (i1==1) then
           if (is_pol==1) call msg('r','[up] '//what//' ',(/(E%f(ib,ik,1),ib=i1,min(i1+7,E%nb))/))
           if (is_pol==2) call msg('r','[dn] '//what//' ',(/(E%f(ib,ik,2),ib=i1,min(i1+7,E%nb))/))
         else
           if (is_pol==1) call msg('r','     '//what//' ',(/(E%f(ib,ik,1),ib=i1,min(i1+7,E%nb))/))
           if (is_pol==2) call msg('r','     '//what//' ',(/(E%f(ib,ik,2),ib=i1,min(i1+7,E%nb))/))
         endif
       endif
     enddo
   enddo
   !
 enddo
 !
end subroutine
