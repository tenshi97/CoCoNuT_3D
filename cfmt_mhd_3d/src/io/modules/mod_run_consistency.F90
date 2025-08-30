module run_consistency
  use precision

  implicit none
! LOCAL variables that are not in modules

  private
  public check_run_consistency

contains

  subroutine check_run_consistency(calculation)
    use precision
    use abort
    use mpi_vertex, only : myproc
    use eos_sn2, only : lsrolo
    use param_rt

    use configure
    use print_stdout_mod
    implicit none
    character*(*),intent(in) :: calculation

    if (trim(calculation) .eq. 'override')    then 
       call printit_taskX(0,"Manual override of consistency check.")
    else if (trim(calculation) .eq. 'hydro')       then
       if (config%p_ntr .ne. 0) raise_abort("check_run_consistency(): hydro mode: wrong config%p_ntr")
       if (config%p_nbk .ne. 0) raise_abort("check_run_consistency(): hydro mode: wrong p_nbk")
    else if (trim(calculation) .eq. 'collapse')    then 
       config%max_central_density = 1.4e14_rk
       if (config%laghyd .ne. 1) raise_abort("check_run_consistency(): collapse mode: wrong laghyd")
       if (config%p_ntr .ne. 1) raise_abort("check_run_consistency(): collapse mode: wrong config%p_ntr")
       if (config%p_nbk .ne. 1) raise_abort("check_run_consistency(): collapse mode: wrong config%p_nbk")
       !if (rolo_den .ne. 1._rk) raise_abort("check_run_consistency(): collapse mode: wrong rolo_den")
       if (config%eos_sw .ne. 0) raise_abort("check_run_consistency(): collapse mode: wrong eos_sw")

       if (lsrolo .gt. 1.e10_rk) raise_abort("check_run_consistency(): collapse mode: too large lsrolo")
    else if (trim(calculation) .eq. 'bounce')      then
       config%stop_bounce = 1
       if (config%laghyd .ne. 0) raise_abort("check_run_consistency(): bounce mode: wrong laghyd")
       if (config%p_ntr .ne. 1) raise_abort("check_run_consistency(): bounce mode: wrong config%p_ntr")
       if (config%p_nbk .ne. 1) raise_abort("check_run_consistency(): bounce mode: wrong config%p_nbk")
       !if (rolo_den .ne. 1._rk) raise_abort("check_run_consistency(): bounce mode: wrong rolo_den")
       if (config%eos_sw .ne. 0) raise_abort("check_run_consistency(): bounce mode: wrong eos_sw")

       if (lsrolo .gt. 1.e10_rk) raise_abort("check_run_consistency(): bounce mode: too large lsrolo")
       
    else if (trim(calculation) .eq. 'postbounce')  then 
       
       if (config%laghyd .ne. 0) raise_abort("check_run_consistency(): postbounce mode: wrong laghyd")
       if (config%p_ntr .ne. 1) raise_abort("check_run_consistency(): postbounce mode: wrong config%p_ntr")
       if (config%p_nbk .ne. 1) raise_abort("check_run_consistency(): postbounce mode: wrong config%p_nbk")
    else
       raise_abort("check_run_consistency(): calculation mode unknown")
    end if

  end subroutine check_run_consistency
  
end module run_consistency
