
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
! Partialy inspired from Octopus v. 3.1.0
!
! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
!
subroutine NL_linear_response(en)
 !
 use YPP_non_linear, ONLY:NL_response,Use_FFT,l_force_min_damp,W_minimal_step
 use electrons,      ONLY:levels
 use pars,           ONLY:SP,DP,cI,cZERO
 use nl_optics,      ONLY:NL_P_t,NL_step
 use units,          ONLY:pi,HA2EV
 use electrons,      ONLY:levels
 use vec_operate,    ONLY:v_norm
 use memory_m,       ONLY:mem_est
 use fields,         ONLY:Efield,small_a_frequency,small_a
 use functions,      ONLY:theta_function
 use units,          ONLY:AUT2FS
 use com,            ONLY:error,msg
 !
 implicit none
 !
 type(levels)  :: en
 !
 ! Work Space
 !
 integer  :: i_t,i_w,i_p_dir,iE_err,i_Probe,i_field,n_fields
 character(10) :: FT_title
 real(SP)      :: W_step,W_reached,Field_time,outdata(7),E_ver(3,3)
 complex(SP)   :: Pol_proj(3),E_w_SP(3),Field,E_of_t(2),A_of_t(2),theta,delta,Field_freq
 complex(SP), allocatable :: EELS(:,:),EPS(:,:),freqs(:)
 complex(DP), allocatable :: Pol_w(:,:),Pol_t(:,:),E_w(:,:),E_t(:,:)
 !
 ! call section("=",'Optics')
 !=========================
 !
 i_Probe=1
 n_fields=1
 !
 ! FFT setup
 !--------------
 !
 allocate(EELS(NL_response%W_n_steps,3),EPS(NL_response%W_n_steps,3),Pol_w(NL_response%W_n_steps,3)) 
 allocate(E_w(NL_response%W_n_steps,n_fields),Pol_t(NL_response%T_n_steps,3))
 !
 call mem_est("EELS EPS Pol_w",(/size(EELS),size(EPS),size(Pol_w)/),(/2*SP,2*SP,2*DP/))
 call mem_est("E_w Pol_t",(/size(E_w),size(Pol_t)/),(/2*DP,2*DP/))
 !
 EELS   =cZERO
 EPS    =cZERO
 Pol_w  =cZERO
 Pol_t  =cZERO
 E_w    =cZERO
 !
 ! Fourier Transform 
 !====================
 !
 W_step=(NL_response%W_range(2) - NL_response%W_range(1))/(NL_response%W_n_steps-1._SP)
 !
 allocate(freqs(NL_response%W_n_steps))
 freqs=cZERO
 !
 do i_w = 1, NL_response%W_n_steps
   freqs(i_w)=dcmplx((i_w-1)*W_step+NL_response%W_range(1),0._SP)
 enddo
 !
 ! Electric Field...
 !-------------------
 !
 iE_err=0
 do i_w = 1, NL_response%W_n_steps
   !
   if(trim(NL_response%damping) /= "NONE")                      Field_freq=freqs(i_w)-cI*NL_response%damp_factor
   if(trim(NL_response%damping) == "NONE".and.l_force_min_damp) Field_freq=freqs(i_w)-cI*NL_response%Field_damp_factor
   !
   E_w(i_w,1)=Efield(i_Probe)%amplitude*small_a_frequency(Field_freq,Efield(i_Probe),iE_err)
   E_w(i_w,1)=E_w(i_w,1)*exp(cI*freqs(i_w)*NL_step)  
   !
   if (iE_err/=0) call error("Error in small_a_frequency!")
   !
 enddo
 !
 !... FT
 !
 allocate(E_t(NL_response%T_n_steps,n_fields))
 do i_field=1,n_fields
   do i_t=1,NL_response%T_n_steps
     Field_time=(i_t-Efield(i_Probe)%t_initial_indx)*NL_step
     A_of_t =small_a(Field_time,Efield(i_Probe),0)
     E_of_t =small_a(Field_time,Efield(i_Probe),1)
     theta  =theta_function(Field_time,NL_step,0)
     delta  =theta_function(Field_time,NL_step,1)
     E_t(i_t,i_field)=Efield(i_Probe)%amplitude*(E_of_t(i_field)*theta+A_of_t(i_field)*delta)
   enddo
   !
   if(trim(NL_response%damping) == "NONE".and.l_force_min_damp) then
     call NL_Cdamp_it("LORENTZIAN",NL_response%Field_damp_factor,E_t(:,i_field),1)
   else
     call NL_Cdamp_it(NL_response%damping,NL_response%damp_factor,E_t(:,i_field),1)
   endif
   !
   FT_title=' E#1'
   !
   call NL_1D_Fourier_Transform('T2W',FT_title,freqs,E_w(:,i_field),NL_response%W_n_steps,&
&                                E_t(:,i_field),NL_response%T_n_steps,1)
   !
 enddo
 deallocate(E_t)
 !
 do i_p_dir=1,3
   !
   if (i_p_dir==1) FT_title=' P (X axis)'
   if (i_p_dir==2) FT_title=' P (Y axis)'
   if (i_p_dir==3) FT_title=' P (Z axis)'
   !
   ! Polarization
   !  
   Pol_t(:,i_p_dir)=NL_P_t(NL_response%i_T_range(1):NL_response%i_T_range(2),i_p_dir)
   call NL_1D_Fourier_Transform('T2W',FT_title,freqs,Pol_w(:,i_p_dir),NL_response%W_n_steps,&
&                               Pol_t(:,i_p_dir),NL_response%T_n_steps,1)
  !
 enddo
 !
 ! a) If the induced field is     computed: E_input=E_ext
 ! b) If the induced field is NOT computed: E_input=E_tot
 !
 ! a) P(w,ic) = E_ext(w,ic') alpha_RID(w,ic',ic) = |E_ext(w)| E_versor(ic') alpha_RID(w,ic,ic')
 ! b) P(w,ic) = E_tot(w,ic') alpha_IRR(w,ic',ic) = |E_tot(w)| E_versor(ic') alpha_IRR(w,ic,ic')
 ! 
 ! a) alpha_RID(w,ic,ic') = ( P(w,ic),E_versor(ic') )/|E_ext(w)| 
 ! b) alpha_IRR(w,ic,ic') = ( P(w,ic),E_versor(ic') )/|E_tot(w)|
 !
 ! Then I use 
 ! a) EPS_{ic,ic'}(w) = (1-4 pi alpha_RID(w,ic))^{-1}
 ! b) EPS_{ic,ic'}(w) =  1+4 pi alpha_IRR(w,ic)
 !
 ! to calculate EELS I use
 !
 ! EPS=-1./EELS
 !
 if(Efield(i_Probe)%ef_pol=="linear"  ) then
   !      
   E_ver(:,1) = Efield(i_Probe)%versor
   E_ver(:,2)=(/0,1,0/)
   E_ver(:,3)=(/0,0,1/)
   !
   if( (1._SP-dot_product( E_ver(:,1),E_ver(:,2) ))<1.E-5) E_ver(:,2)=(/1,0,0/)
   E_ver(:,2) = E_ver(:,2)-dot_product( E_ver(:,1),E_ver(:,2) )*E_ver(:,1)
   E_ver(:,2) = E_ver(:,2)/v_norm(E_ver(:,2))
   !
   if( (1._SP-dot_product( E_ver(:,1),E_ver(:,3) ))<1.E-5) E_ver(:,3)=(/1,0,0/)
   E_ver(:,3) = E_ver(:,3)-dot_product( E_ver(:,1),E_ver(:,3) )*E_ver(:,1) &
&                         -dot_product( E_ver(:,2),E_ver(:,3) )*E_ver(:,2)
   E_ver(:,3) = E_ver(:,3)/v_norm(E_ver(:,3))
   !
 endif
 !
 do i_w = 1, NL_response%W_n_steps
   !
   if(Efield(i_Probe)%ef_pol=="linear"  ) Field = cmplx(E_w(i_w,1),kind=SP)
   !
   do i_p_dir=1,3
     !
     Pol_proj(i_p_dir)=dot_product(E_ver(:,i_p_dir),cmplx(Pol_w(i_w,:),kind=SP) )
     EPS( i_w,i_p_dir) = 4._SP*pi*Pol_proj(i_p_dir)/Field
     !
   enddo
   !
 enddo
 !
 EPS(:,1)  =  1._SP+EPS(:,1)
 EELS(:,:) = -1._SP/EPS(:,:)
 !
 ! Write Output
 !
 call LR_write_output()
 !
 deallocate(EELS,EPS,Pol_w,Pol_t,E_w)
 call mem_est("EELS EPS Pol_w Pol_t E_w")
 !
 return
 !
 contains
   !      
   subroutine LR_write_output()
     !
     use com,            ONLY:msg,of_open_close
     use stderr,         ONLY:set_real_printed_length 
     use pars,           ONLY:schlen
     !
     implicit none
     !
     ! Work Space
     !
     character(schlen) :: eps_file_name,eels_file_name,efield_file_name, &
&                         pol_file_name
     character(schlen) :: headings(11)
     integer           :: i_t
     !
     call set_real_printed_length(f_length=15,g_length=15)
     !
     eps_file_name  ="YPP-eps_along_E"
     eels_file_name ="YPP-eels_along_E"
     efield_file_name="YPP-efield_along_E"
     !
     headings(1)="E [eV]"
     !
     call of_open_close(eps_file_name,'ot')
     headings(2:7)=(/'Im/eps_d1','Re/eps_d1','Im/eps_d2','Re/eps_d2','Im/eps_d3','Re/eps_d3'/)
     call msg('o '//trim(eps_file_name),'#')
     call msg('o '//trim(eps_file_name),'#',headings(1:7),INDENT=0,USE_TABS=.TRUE.)
     call msg('o '//trim(eps_file_name),'#')
     !
     call of_open_close(eels_file_name,'ot') 
     headings(2:7)=(/'Im/eel_d1','Re/eel_d1','Im/eel_d2','Re/eel_d2','Im/eel_d3','Re/eel_d3'/)
     call msg('o '//trim(eels_file_name),'#')
     call msg('o '//trim(eels_file_name),'#',headings(1:7),INDENT=0,USE_TABS=.TRUE.)
     call msg('o '//trim(eels_file_name),'#')
     !
     call of_open_close(efield_file_name,'ot')
     headings(2:7)=(/"E_x/Re[V/m]","E_x/Im[V/m]","E_y/Re[V/m]","E_y/Im[V/m]","E_z/Re[V/m]","E_z/Im[V/m]"/)
     call msg('o '//trim(efield_file_name),'#')
     call msg('o '//trim(efield_file_name),'#',headings(1:7),INDENT=0,USE_TABS=.TRUE.)
     call msg('o '//trim(efield_file_name),'#')
     !
     do i_w = 1, NL_response%W_n_steps
      !
      W_reached = real(freqs(i_w))*HA2EV
      outdata(1)=W_reached
      !
      E_w_SP = cZERO
      if(Efield(i_Probe)%ef_pol=="linear"  ) E_w_SP=Efield(i_Probe)%versor*cmplx(E_w(i_w,1),kind=SP)
      !
      outdata(2:4)=(/ real(E_w_SP(1) ),aimag(E_w_SP(1) ),real(E_w_SP(2) ) /)
      outdata(5:7)=(/ aimag(E_w_SP(2)),real(E_w_SP(3)),aimag(E_w_SP(3)) /)
      call msg('o '//trim(efield_file_name),'',outdata(1:7),USE_TABS=.TRUE.)
      !
      outdata(2:4)=(/ aimag(EPS(i_w,1) ),real(EPS(i_w,1) ),aimag(EPS(i_w,2) ) /)
      outdata(5:7)=(/ real( EPS(i_w,2) ),aimag( EPS(i_w,3) ),real( EPS(i_w,3) ) /)
      call msg('o '//trim(eps_file_name),'',outdata(1:7),USE_TABS=.TRUE.)
      !
      outdata(2:4)=(/ aimag(EELS(i_w,1) ), real(EELS(i_w,1) ),aimag(EELS(i_w,2) ) /)
      outdata(5:7)=(/ real( EELS(i_w,2) ),aimag(EELS(i_w,3) ),real( EELS(i_w,3) ) /)
      call msg('o '//trim(eels_file_name),'',outdata(1:7),USE_TABS=.TRUE.)
      !
    enddo
    !
    if(trim(NL_response%damping) /= "NONE") then
       !     
       pol_file_name   ="YPP-damped_polarization"
       !
       headings(1)="TimeE [fs]"
       !
       call of_open_close(pol_file_name,'ot')
       headings(2:4)=(/'Pol_x','Pol_y','Pol_z'/)
       call msg('o '//trim(pol_file_name),'#')
       call msg('o '//trim(pol_file_name),'#',headings(1:4),INDENT=0,USE_TABS=.TRUE.)
       call msg('o '//trim(pol_file_name),'#')
       !
       do i_t = NL_response%i_T_range(1),NL_response%i_T_range(2)
         outdata(1)=(i_t-1.0)*NL_step*AUT2FS
         outdata(2:4)=(/ NL_P_t(i_t,1), NL_P_t(i_t,2),NL_P_t(i_t,3) /)
         call msg('o '//trim(pol_file_name),'',outdata(1:7),USE_TABS=.TRUE.)
       enddo
       !
       call of_open_close(pol_file_name)
       !
    endif
    !
    call of_open_close(eps_file_name)
    call of_open_close(eels_file_name)
    call of_open_close(efield_file_name)
    call set_real_printed_length()
   end subroutine LR_write_output
  !
end subroutine NL_linear_response
