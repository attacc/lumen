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
module LOGO
 ! 
 implicit none
 !
 integer      :: code_version(3)
 integer      :: code_GPL_revision
 integer      :: code_revision
 character(50):: code_hash
 !
 integer                     :: n_logo_lines
 integer                     :: ID_logo
 integer                     :: ID_logo_stderr
 integer, parameter, private :: max_n_logo_lines=100
 character(60)               :: logo_line(max_n_logo_lines)
 !
 contains
   !
   integer function pickup_a_random(rand_range)
     use pars,  ONLY:DP,SP
     real(SP) :: rand_range
     !
     ! Work Space
     !
     character(12)      :: ch(3)
     integer            :: iseed(8)
     real(DP), external :: dlaran
     call date_and_time(ch(1),ch(2),ch(3),iseed)
     iseed=iabs(iseed)
     !
     ! iseed(4) must be odd
     !
     iseed(4)=2*(iseed(4)/2)+1
     !
     pickup_a_random=int(dlaran(iseed(4:))*rand_range)
   end function
   !
   subroutine pickup_a_logo(unit_)
     !
     use pars,     ONLY:schlen,SP
     use stderr,   ONLY:intc,slash
     !
     integer      :: unit_
     ! 
     ! Work Space
     !
     character(schlen) :: comment(3),ch,MPI_string,OMP_string
     integer           :: max_length,i1,left_space,iA,ic
     integer           :: n_logos,n_logos_stderr
     !
     ! Get version and revision
     !
     include 'version.inc'
     !
     !
     comment(1)='GPL Version '//&
&                trim(intc(code_version(1)))//'.'//trim(intc(code_version(2)))//&
&                '.'//trim(intc(code_version(3)))//' Revision '//trim(intc(code_GPL_revision))
     !
     MPI_string="Serial"
#if defined _MPI
     MPI_string="MPI"
#endif
     OMP_string=" "
#if defined _OPENMP
     OMP_string="+OpenMP"
#endif
     comment(2)=trim(MPI_string)//trim(OMP_string)//" Build"
     !
     comment(3)='(based on Yambo 4.1.2)'
     !
     n_logos       =1
     n_logos_stderr=1
     !
     if (ID_logo<0) then
       ID_logo=pickup_a_random(real(n_logos,SP))
       ID_logo_stderr=n_logos+pickup_a_random(real(n_logos_stderr,SP))
       if (ID_logo==0) ID_logo=1
       if (ID_logo_stderr==0) ID_logo_stderr=1
     endif
     !
     if (unit_==6) i1=ID_logo_stderr
     if (unit_/=6) i1=ID_logo
     !
     n_logo_lines=9
     !
     select case(i1)
       case(1)
         logo_line( 1)='                                                      '
         logo_line( 2)='//$$        /$$       /$$       $$ /$$$$$$$   $$    $$'
         logo_line( 3)='| $$       //$$    $$| $$$    /$$$| $$       /$$$$  $$'
         logo_line( 4)='| $$       | $$    $$| $$$$  /$$$$| $$      | $$$$$ $$'
         logo_line( 5)='| $$       | $$    $$| $$ $$/$$ $$| $$$$$$$ | $$ $$$$$'
         logo_line( 6)='| $$       | $$    $$| $$  $$$| $$| $$      | $$   $$$'
         logo_line( 7)='| $$       | $$    $$| $$S  $ | $$| $$      | $$    $$'
         logo_line( 8)='| $$$$$$$  | $$$$$$$$| $$     | $$| $$$$$$$/| $$    $$'
         logo_line( 9)='|______/   |________/|__/     |__/|_______/  S__/  S_/'
       case(2)
         logo_line( 1)='                                             '
         logo_line( 2)=' __      _     _  __   __  _______         _ '
         logo_line( 3)='|  |    | |   | ||  |_|  ||  _____||\\     | |'
         logo_line( 4)='|  |    | |   | ||       || |      | \\    | |'
         logo_line( 5)='|  |    | |   | ||       || |____  |\\ \\   | |'
         logo_line( 6)='|  |    | |   | ||       ||  ____| | \\ \\  | |'
         logo_line( 7)='|  |__  | |___| || ||_|| || |_____ | |\\ \\ | |'
         logo_line( 8)='|_____| |_______||_|   |_||_______|| | \\____|'
         logo_line( 9)='                                             '
     end select
     !
     ! Replace "S" with "\" and find the max length of 
     ! the lines to center the two comments
     !
     do i1=1,n_logo_lines
       iA=index(logo_line(i1),'S')
       if (iA==0) cycle
       ch=logo_line(i1)
       do while (index(ch,'S')/=0) 
         iA=index(ch,'S')
         ch(iA:iA)=slash
       enddo
       logo_line(i1)=ch
     enddo
     !
     max_length=-1
     do i1=1,n_logo_lines
       max_length=max( len_trim(logo_line(i1)) , max_length)
     enddo
     !
     n_logo_lines=n_logo_lines+1
     logo_line(n_logo_lines)=' '
     !
     do ic=1,3
       n_logo_lines=n_logo_lines+1
       logo_line(n_logo_lines)=trim(comment(ic))
       left_space=(max_length-len_trim(comment(ic)))/2+1
       if (left_space>0) then
         write (logo_line(n_logo_lines),'(t'//trim(intc(left_space))//',a)') trim(comment(ic))
       endif
     enddo
     !
   end subroutine
   !
end module LOGO
