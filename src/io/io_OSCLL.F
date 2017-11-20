!
!        Copyright (C) 2000-2015 the YAMBO team
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
integer function io_OSCLL(q,ID)
 !
 use pars,          ONLY:SP,schlen
 use X_m,           ONLY:X_t,X_alloc
 use electrons,     ONLY:levels
 use matrix_operate,ONLY:mat_c2r,mat_r2c
 use R_lattice,     ONLY:nXkibz,bz_samp
 use QP_m,          ONLY:QP_nk
 use IO_m,          ONLY:io_connect,io_disconnect,io_sec,frag_DIP,&
&                        io_bulk,read_is_on,write_is_on,io_header,io_extension,&
&                        RD_CL_IF_END,RD_CL,io_action,&
&                        IO_NOT_ALLOWED,io_DIP,io_fragment
 use IO_m,          ONLY:io_elemental,io_status
 use nl_optics,     ONLY:NL_bands,OSCLL,OSCLL_k,ng_oscll,NL_nbands
 use memory_m,      ONLY:mem_est
 use parallel_m,    ONLY:PAR_Xk_ibz_index,PAR_Xk_nibz,myid
#if defined _TIMING
 use timing_m,      ONLY:timing
#endif
 !
 implicit none
 integer     ::ID
 type(bz_samp), intent(in) :: q
 !
 ! Work Space
 !
 integer :: ik,ik_mem,i_fragment,i_sp_pol,iqbz,ig,ib
 character(schlen)    :: VAR_name
 real(SP),allocatable :: disk_data(:,:,:)
 !
#if defined _TIMING
 call timing('io_OSCLL',OPR='start')
#endif
 !
 io_extension(ID)='OSCLL'
 !
 io_OSCLL=io_connect(desc=trim(io_extension(ID)),type=2,ID=ID)
 if (io_OSCLL/=0) goto 1
 !
 if (any((/io_sec(ID,:)==1/))) then
   !
   io_OSCLL=io_header(ID,R_LATT=.true.,WF=.true.,IMPOSE_SN=.true.,T_EL=.true.)
   if (io_OSCLL/=0) goto 1
   !
   ! Variable for the Oscillators
   !
   call io_elemental(ID,VAR="OSCLLPARS",VAR_SZ=3,MENU=0)
   call io_elemental(ID,VAR="[NL] RL vectors in oscillators :",I0=ng_oscll,CHECK=.true.,OP=(/"<="/))
   call io_elemental(ID,VAR="[NL] OSCLL band rane           :",I1=NL_bands,CHECK=.true.,OP=(/"==","=="/))
   call io_elemental(ID,VAR="",VAR_SZ=0)
   !
   io_OSCLL=io_status(ID)
   if (io_OSCLL/=0) goto 1
   !
 endif
 !
 i_fragment=maxval(io_sec(ID,:))-1
 !
 if(i_fragment<=nXkibz) then
   i_sp_pol=1
   ik=i_fragment
 else
   i_sp_pol=2
   ik=i_fragment-nXkibz
 endif
 !
 if (ik==0) goto 1
 !
 ik_mem=PAR_Xk_ibz_index(ik)
 !
 allocate(disk_data(NL_nbands,q%nbz,2))
 call mem_est("disk_data",(/size(disk_data)/),(/SP/))
 !
 if (read_is_on(ID)) then
   if(.not.allocated(OSCLL)) then
     allocate(OSCLL(NL_bands(1):NL_bands(2),NL_bands(1):NL_bands(2),PAR_Xk_nibz,q%nbz,ng_oscll))
     call mem_est("OSCLL",(/size(OSCLL)/),(/2*SP/))     
   endif
 endif
 !
 ! Fragmentation
 !
 call io_fragment(ID,i_fragment=i_fragment)
 !
 ! Manage RD_CL_IF_END
 !
 if (io_action(ID)==RD_CL_IF_END.and.ik==nXkibz) io_action(ID)=RD_CL
 !
 do ib=NL_bands(1),NL_bands(2)
   do ig=1,ng_oscll
     !
     if (write_is_on(ID)) call mat_c2r(OSCLL_k(ib,:,:,ig),disk_data)
     !
     write (VAR_name,'(3(a,i4.4))') 'OSCLL_',ik,'_ib_',ib,'_ig_',ig
     call io_bulk(ID,VAR=trim(VAR_name),VAR_SZ=shape(disk_data))
     call io_bulk(ID,R3=disk_data)
     !
     io_OSCLL=io_status(ID)
     if (io_OSCLL/=0) goto 1
     !
     if (read_is_on(ID)) call mat_r2c(disk_data,OSCLL(ib,:,ik_mem,:,ig))
     !
   enddo
 enddo
 !
 deallocate(disk_data)
 call mem_est("disk_data")
 !
1 call io_disconnect(ID=ID)
 !
#if defined _TIMING
 call timing('io_OSCLL',OPR='stop')
#endif
 !
end function
