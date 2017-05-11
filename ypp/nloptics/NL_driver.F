!
!        Copyright (C) 2000-2015 the YAMBO team
!              http://www.yambo-code.org
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
subroutine NL_driver(en)
 !
 use YPP_non_linear, ONLY:X_order,NL_response,W_minimal_step
 use pars,           ONLY:SP,pi,cZERO,lchlen
 use units,          ONLY:HA2EV,AUT2FS
 use electrons,      ONLY:levels
 use parser_m,       ONLY:parser
 use nl_optics,      ONLY:n_runs,NL_P_t,NL_damping,NL_steps,NL_step, &
&                         NL_steps,n_intensities
 use fields,         ONLY:Efield,EtoT,Divide_by_Field,ext_field,Efield_duplicate,Efield_strength
 use IO_m,           ONLY:OP_RD_CL,DUMP,io_control,NONE
 use com,            ONLY:msg,error,warning
 use stderr,         ONLY:intc,real2ch
 use memory_m,       ONLY:mem_est
 use LIVE_t,         ONLY:live_timing
 !
 implicit none
 !
 type(levels)  :: en
 !
 ! Working Space
 !
 integer, external::io_NL
 integer  :: io_NonLinear,ID,i_Probe
 integer  :: i_r,i_d,i_order,i_t_field_range(2),i_time_range_save(2)
 real(SP) :: W_step,T_period,time_range_save(2)
 real(SP),    allocatable  :: NL_P(:,:,:),Harmonic_Frequency(:,:),Probe_frequency(:)
 complex(SP), allocatable  :: X_effective(:,:,:),Susceptibility(:,:,:)
 logical  :: Probe_consistency_error
 type(ext_field)   :: Efield_probe
 real(SP) :: NL_tot_time
 !
 call section("*",'NL Post Processing')
 !
 ! Input/Output
 !
 call io_control(ACTION=OP_RD_CL,COM=NONE,SEC=(/1/),MODE=DUMP,ID=ID)
 io_Nonlinear=io_NL(ID)
 if(io_Nonlinear/=0) call error("Error reading NonLinear database")
 !
 allocate(NL_P(NL_steps,n_runs,3),Probe_frequency(n_runs))
 !
 ! Reads the polarization from dbs containing simulation at different frequencies 
 !
 Probe_consistency_error=.FALSE.
 i_Probe=1
 !
 do i_r=1,n_runs
   call io_control(ACTION=OP_RD_CL,COM=NONE,SEC=(/1+i_r/),MODE=DUMP,ID=ID)
   io_Nonlinear=io_NL(ID)
   !
   Efield(i_Probe)%amplitude=Efield_strength(Efield(i_Probe)%intensity,'AU')
   !
   if (i_r==1) then
     call Efield_duplicate(Efield(i_Probe),Efield_probe)
   else
     if (Efield_probe%t_initial/=Efield(i_Probe)%t_initial)   Probe_consistency_error=.TRUE.
     if (Efield_probe%ef_name  /=Efield(i_Probe)%ef_name)     Probe_consistency_error=.TRUE.
   endif
   !
   do i_d=1,3
     NL_P(:,i_r,i_d)=NL_P_t(:,i_d)
   enddo
   !
   Probe_frequency(i_r)=Efield(i_Probe)%frequency(1)
   !
 enddo
 !
 Efield_probe%ef_name=trim(Efield_probe%ef_name)
 !
 if (Probe_consistency_error) call error('Not all probes have consistent parameters')
 !
 NL_tot_time    =real(NL_steps,SP)*NL_step
 !
 ! Integer time range
 !
 if(NL_response%T_range(1)>0._SP) then
   NL_response%i_T_range(1)=nint(NL_response%T_range(1)/NL_step)+1
 else
   NL_response%i_T_range(1)=1
   NL_response%T_range(1)  =0._SP
 endif
 !
 if(NL_response%T_range(2)>0._SP) then 
   NL_response%i_T_range(1)=nint(NL_response%T_range(2)/NL_step)
 else
   NL_response%i_T_range(2)=NL_steps
   NL_response%T_range(2)  =(NL_steps-1)*NL_step
 endif
 !
 if(Efield_probe%ef_name=='SIN'.or.Efield_probe%ef_name=='SOFTSIN'.or.Efield_probe%ef_name=='ANTIRES') then
   !        
   ! Time Initialization for periodic probes
   !
   W_step         =minval(Probe_frequency) ! smaller frequency
   T_period       =EtoT(W_step)            ! largest period
   !
   ! Check if the time for dephasing is enought
   ! 
   if(NL_response%T_range(1)<=2._SP/NL_damping*5._SP) &
&     call warning('Dephasing time too short, set initial time >= '//trim(real2ch(2._SP/NL_damping*5._SP*AUT2FS))//' fs')
   !
   ! Check if total simulation time is long enough for analysis  
   !
   if((NL_response%T_range(2) - NL_response%T_range(1))<T_period) &
&     call warning(' Simulation time not long enough. Period is'//trim(real2ch(T_period*AUT2FS))//' fs')
   !
 endif
 !
 call msg('s','[NL] Field type    :', trim(Efield_probe%ef_name))
 call msg('s','[NL] Integer time range :',NL_response%i_T_range)
 call msg('s','[NL] Time range :',NL_response%T_range*AUT2FS)
 !
 NL_response%T_n_steps=NL_response%i_T_range(2)-NL_response%i_T_range(1)+1
 !
 ! Defined the smallest frequency step from the simulation length
 !
 W_minimal_step=2._SP*pi/(real(NL_response%T_n_steps-1._SP,SP))*NL_step
 ! 
 select case (trim(Efield_probe%ef_name))
   !
   ! CONTINUOUS MONOCHROMATIC SOURCES
   !
 case('SIN','SOFTSIN','ANTIRES')
     !  
     allocate(X_effective(0:X_order,n_runs,3))
     allocate(Susceptibility(0:X_order,n_runs,3),Harmonic_Frequency(0:X_order,n_runs))
     !
     do i_order=0,X_order
       Harmonic_Frequency(i_order,:)=real(i_order,SP)*Probe_Frequency(:)
     enddo
     !
     call msg('s',':: Harmonic frequencies    [eV]:',(/minval(Harmonic_Frequency),maxval(Harmonic_Frequency)/)*HA2EV)
     !
     ! Find the Fourier coefficients by inversion
     !
     call msg('s','[NL] Effective max time period [fs]:',T_period*AUT2FS)
     !
     call live_timing('X effective :', n_runs)
     !
     time_range_save  =NL_response%T_range
     i_time_range_save=NL_response%i_T_range
     !
     do i_r = 1, n_runs
       !
       T_period=T_range(1,i_r)
       !
       do i_d = 1,3
         call NL_coefficients_Inversion(X_order+1,X_order+1,NL_P(:,i_r,i_d),X_effective(0:,i_r,i_d),&
&                                     Harmonic_Frequency(0:,i_r),T_period)
       end do
       !
       call live_timing(steps=1)
       !
     enddo
     !
     call live_timing()
     !
     NL_response%T_range=time_range_save
     NL_response%i_T_range=i_time_range_save
     !
     ! Calculate Susceptibilities from X_effective
     !
     do i_order=0,X_order
       !
       do i_r=1,n_runs
         !
         if (i_order==1) then
           Susceptibility(i_order,i_r,1)=4._SP*pi*dot_product( Efield(i_Probe)%versor(:),X_effective(i_order,i_r,:) )*&
                &                                                   Divide_by_Field(Efield(i_Probe),i_order)
         else 
           Susceptibility(i_order,i_r,:)=4._SP*pi*X_effective(i_order,i_r,:)*Divide_by_Field(Efield(i_Probe),i_order)
         end if
         !
       enddo
       ! 
     enddo
     !
     ! Richardson extrapolation
     !
     if(n_intensities>1) then
        call NL_richardson(Susceptibility)
     endif
     !
     ! Plot
     !
     call plot_X_effective()
     deallocate(X_effective,Susceptibility,Harmonic_Frequency)
     !
  case ('DELTA')
    !
    ! For testing purpose
    !
    if(trim(NL_response%damping) /= "NONE") then
       call NL_Rdamp_it(NL_response%damping,NL_response%damp_factor,&
&                NL_P_t(NL_response%i_T_range(1):NL_response%i_T_range(2),1:3),3)
    endif
    !
    call NL_linear_response(en)
    deallocate(NL_P_t)
     !
     ! OTHER CASES NOT YET WORKING
     !
  case default 
    call error('Analysis not available for the '//trim((Efield_probe%ef_name))//' field')
  end select
  !
  deallocate(NL_P)
  !
contains
   !
   real(SP) function T_range(i_n,i_f)
     !
     integer  :: i_n,i_f
     real(SP) :: dephasing_time
     !
     dephasing_time=NL_response%T_range(1)
     !
     NL_response%T_range=time_range_save
     !
     if (i_f>  0) T_range = EtoT(E=Harmonic_Frequency(i_n,i_f))
     if (i_f== 0) T_range = EtoT(E=maxval(Harmonic_Frequency(:,:)))
     !
     if (i_f==-1) then
       T_range = NL_response%T_range(2)- NL_response%T_range(1)
       return
     endif
     !
     NL_response%T_range(2)  = NL_steps*NL_step
     NL_response%T_range(1)  = NL_response%T_range(2) - T_range
     !
     if (NL_response%T_range(1)< dephasing_time) then
       call warning(" For frequency number "//intc(i_f)//" simulation length too short ")
     endif 
     !
     NL_response%i_T_range(1)=nint(NL_response%T_range(1)/NL_step)+1
     NL_response%i_T_range(2)=nint(NL_response%T_range(2)/NL_step)
     !
   end function 
   !
end subroutine
