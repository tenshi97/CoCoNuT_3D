program dump_test
  use testdumper
  use precision
  real(kind=rk) :: foo(3,4)
  real(kind=rk) :: dice

  foo = 1.0_rk

  call random_seed()
  call random_number(dice)

  if (dice >= 0.5) then
    foo(2,3) = dice
  endif

  write(*,*) "dice = ", dice

  dump(foo)

end program dump_test
