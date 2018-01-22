!
!        Copyright (C) 2000-2017 the LUMEN team
!
! Authors (see AUTHORS file for details): MG CA
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
subroutine NL_test_collisions(k,q) 
 !
 ! Test-collisions built on the fly vs collisions stored in memory
 !
 use pars,           ONLY:cZERO,SP
 use R_lattice,      ONLY:bz_samp,qindx_S
 use collision,      ONLY:Fock_collisions,SEX_collisions,nFock_collisions,nSEX_collisions
 use QP_m,           ONLY:QP_table,QP_nk,QP_n_states
 use nl_optics,      ONLY:NL_bands
 !
 implicit none
 ! 
 type(bz_samp), intent(in)  :: k,q
 !
 integer :: i_kmq,i_qp,ib,ibp
 integer :: ic1,ic2
 complex(SP) :: COLL(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),QP_nk)
 complex(SP) :: COLL_OSCLL(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),QP_nk)
 !
 do i_qp=1,QP_n_states
   !
   COLL=cZERO
   !
   call NL_SEX_collisions(k,q,i_qp,COLL_OSCLL) 
   !
   ic1=1
   ic2=1
   !
   do i_kmq=1,QP_nk
     !
     do ibp=NL_bands(1),NL_bands(2)
       do ib=NL_bands(1),NL_bands(2)
         if (nFock_collisions>0) then
           if (Fock_collisions(i_qp)%table(ib,ibp,i_kmq)/=0) then
              COLL(ib,ibp,i_kmq)=COLL(ib,ibp,i_kmq)+Fock_collisions(i_qp)%v_c(ic1)
             ic1=ic1+1
           endif
         endif
         if (nSEX_collisions>0) then
           if (SEX_collisions(i_qp)%table(ib,ibp,i_kmq)/=0) then
             COLL(ib,ibp,i_kmq)=COLL(ib,ibp,i_kmq)+SEX_collisions(i_qp)%v_c(ic2)
             ic2=ic2+1
           endif
         endif
         !
         write(*,*) ' Collision indexex ',ib,ibp,i_kmq
         write(*,*) COLL(ib,ibp,i_kmq)
         write(*,*) COLL_OSCLL(ib,ibp,i_kmq)
         !
       enddo
     enddo
     !
   enddo
   !
   stop
   !
 enddo 
 !
end subroutine NL_test_collisions
