!
!        Copyright (C) 2000-2016 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM CA DS
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
integer function io_COLLISIONS(i_qp,ID,COLL)
 !
 use pars,          ONLY:SP,schlen,lchlen
 use R_lattice,     ONLY:RIM_ng,RIM_n_rand_pts
 use collision,     ONLY:V_kind,COLLISIONS_cutoff,COLLISIONS_element,COLLISIONS_IO_status,&
&                        n_COLLISIONS,COLLISIONS_state,COLL_bands
 use IO_m,          ONLY:io_sec,io_connect,io_header,io_status,io_elemental, &
&                        io_bulk,write_is_on,io_disconnect,&
&                        IO_OUTDATED_DB,read_is_on,io_fragment,IO_mode,&
&                        VERIFY,io_extension,io_raw_extension
 use global_XC,     ONLY:SC_Fock,global_XC_string,EXT_NONE,SE_NONE,SC_COHSEX,SE_GW_NEQ,SE_GW_QED
 use com,           ONLY:error
 use SC,            ONLY:SC_X_size,SC_X_bands
 use QP_m,          ONLY:QP_ng_Sx,QP_n_states,SC_band_mixing,QP_QED_ng
 use parallel_m,    ONLY:PAR_IND_QP
 use X_m,           ONLY:n_X_descs,X_description
 use stderr,        ONLY:intc
 use drivers,       ONLY:l_elel_scatt,l_elphoton_scatt
#if defined _TIMING
 use timing_m,      ONLY:timing
#endif
 !
 implicit none
 !
 integer, intent(inout)            ::ID
 integer, intent(in)               ::i_qp
 type(COLLISIONS_element),optional ::COLL
 !
 ! Work Space
 !
 character(schlen)     :: db_name
 real(SP), allocatable :: COLLISIONS_v_disk(:,:)
 integer               :: V_kind_disk,i_desc,i_qp_loop,io_err
 character(lchlen)     :: V_kind_string,X_string
 !
#if defined _TIMING
 call timing('io_COLLISIONS',OPR='start')
#endif
 !
 if(V_kind==SC_Fock) then
   db_name                   = "COLLISIONS_FOCK"
 elseif(V_kind==SC_COHSEX) then
   db_name                   = "COLLISIONS_COHSEX"
 else
   call error("COLLISIONS I/O : Unknown database type ")
 endif
 !
 io_COLLISIONS=io_connect(desc=db_name,type=2,ID=ID)
 if(io_COLLISIONS/=0) goto 1
 !
 if (any(io_sec(ID,:)==1)) then
   !
   io_COLLISIONS=io_header(ID,QPTS=.true.,R_LATT=.true.,WF=.true.,IMPOSE_SN=.true.)
   !
   call io_elemental(ID,VAR="PARAMETERS",VAR_SZ=10,MENU=0)
   !
!   if (ver_is_gt_or_eq(ID,revision=5152)) then
!     call io_elemental(ID,I0=n_COLLISIONS,&
!&         VAR=' Total number of COLLISIONS      :',CHECK=.true.,OP=(/"=="/))
!   else
   n_COLLISIONS=QP_n_states
!   endif
   !
   call io_elemental(ID,I1=COLL_bands,&
&       VAR=' Bands                           :',CHECK=.true.,OP=(/"==","=="/))
   !
   call io_elemental(ID,R0=SC_band_mixing,&
&       VAR=' Bands mixing                    :',CHECK=.true.,OP=(/"=="/))
   call io_elemental(ID,I0=V_kind,DB_I0=V_kind_disk,CHECK=.true.,OP=(/"=="/))
   call io_elemental(ID,VAR="",VAR_SZ=0,MENU=0)
   !
   V_kind_string     =global_XC_string(V_kind,SE_NONE,0,EXT_NONE)
   !
   call io_elemental(ID,VAR='INTERACTION_KIND',CH0="",VAR_SZ=1,MENU=0)
   call io_elemental(ID,CH0=V_kind_string,&
&       VAR=' Interaction Kind                :',CHECK=.true.,OP=(/"=="/))
   call io_elemental(ID,VAR="",VAR_SZ=0,MENU=0)
   !
   if (V_kind_disk==SC_Fock) then
     !
     call io_elemental(ID,VAR="INTERACTION_PARAMETERS",VAR_SZ=4,MENU=0)
     call io_elemental(ID,I0=QP_ng_Sx,&
&         VAR=' Exchange RL vectors             :',CHECK=.true.,OP=(/"=="/))
     call io_elemental(ID,R0=COLLISIONS_cutoff,&
&         VAR=' Exchange Cutoff                 :',CHECK=.true.,OP=(/"=="/))
     !
     call io_elemental(ID,VAR="",VAR_SZ=0)
     !
   else if (V_kind_disk==SC_COHSEX) then
     !
     ! Global Response Descriptions (to be done)
     !==============================
     !
     call io_elemental(ID,VAR="INTERACTION_PARAMETERS",VAR_SZ=2,MENU=0)
     call io_elemental(ID,I0=RIM_ng,VAR=        ' RIM RL components      :',WARN=.true.)
     call io_elemental(ID,I0=RIM_n_rand_pts,VAR=' RIM random points      :',WARN=.true.)
     call io_elemental(ID,VAR="",VAR_SZ=0)
     !
   else
     !
     call io_elemental(ID,VAR="",VAR_SZ=0)
     !
   endif
   !
   ! States
   !
   if (read_is_on(ID).and..not.allocated(COLLISIONS_state)) then
      allocate(COLLISIONS_state(n_COLLISIONS,4))
    endif
   call io_bulk(ID,VAR='COLLISIONS_STATE',VAR_SZ=(/n_COLLISIONS,4/))
   call io_bulk(ID,I2=COLLISIONS_state)
   !
   io_COLLISIONS=io_status(ID)
   if (io_COLLISIONS/=0) goto 1
   !
 endif
 !
 ! Check the restart point
 !
 if (IO_mode(ID)==VERIFY.and.allocated(COLLISIONS_IO_status)) then
   !
   COLLISIONS_IO_status=1
   do i_qp_loop=1,n_COLLISIONS
     if (.not.PAR_IND_QP%element_1D(i_qp_loop)) cycle
     io_extension(ID)=io_raw_extension(ID)
     call io_fragment(ID,i_fragment=i_qp_loop,ierr=io_err)
     if (io_err==0) then
       COLLISIONS_IO_status(i_qp_loop)=1
     else
       COLLISIONS_IO_status(i_qp_loop)=-1
       exit
     endif
   enddo
   !
 endif
 !
 if (any(io_sec(ID,:)==2)) then     
   !
   ! QP-based fragmentation
   !
   call io_fragment(ID,i_fragment=i_qp)
   !
   call io_bulk(ID,VAR='N_COLLISIONS_STATES',VAR_SZ=(/1/))
   call io_bulk(ID,I0=COLL%N)
   !
 endif
 !
 if (any(io_sec(ID,:)==3)) then     
   !
   ! QP-based fragmentation
   !
   if (.not.any(io_sec(ID,:)==2)) call io_fragment(ID,i_fragment=i_qp)
   !
   call io_bulk(ID,VAR='COLLISIONS_TABLE',&
&                    VAR_SZ=(/COLL%D(1,2)-COLL%D(1,1)+1,&
&                             COLL%D(2,2)-COLL%D(2,1)+1,&
&                             COLL%D(3,2)-COLL%D(3,1)+1/))
   !
   call io_bulk(ID,L3=COLL%table(COLL%D(1,1):COLL%D(1,2),&
&                                COLL%D(2,1):COLL%D(2,2),&
&                                COLL%D(3,1):COLL%D(3,2)))
   !
   if(V_kind==SE_GW_NEQ.or.V_kind==SE_GW_QED) then
     ! Real collisions
     call io_bulk(ID,VAR='COLLISIONS_v',VAR_SZ=(/COLL%N/))
     call io_bulk(ID,R1=COLL%v_r)
   else
     ! Complex collisions
     allocate(COLLISIONS_v_disk(COLL%N,2))
     call io_bulk(ID,VAR='COLLISIONS_v',VAR_SZ=(/COLL%N,2/))
     if(write_is_on(ID)) then
       COLLISIONS_v_disk(:,1) =real(COLL%v_c(1:COLL%N))
       COLLISIONS_v_disk(:,2)=aimag(COLL%v_c(1:COLL%N))
       call io_bulk(ID,R2=COLLISIONS_v_disk)
     else
       call io_bulk(ID,R2=COLLISIONS_v_disk)
       COLL%v_c(1:COLL%N)=CMPLX(COLLISIONS_v_disk(:,1),COLLISIONS_v_disk(:,2))
     endif
     deallocate(COLLISIONS_v_disk)
   endif
   !
 endif
 !
1 continue
 !
 call io_disconnect(ID=ID)
 !
#if defined _TIMING
 call timing('io_COLLISIONS',OPR='stop')
#endif
 !
end function io_COLLISIONS
