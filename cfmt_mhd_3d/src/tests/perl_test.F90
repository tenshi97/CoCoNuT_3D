module perl_test_mod
  implicit none

  interface array_info
!PERL-START for(my $dim = 1; $dim < 7; $dim++) {
    module procedure array_info_@[[$dim]]
!PERL-END }
  end interface

  contains

!PERL-START for (my $dim = 1; $dim < 7; $dim++) {
  subroutine array_info_@[[$dim]](a)
    real, dimension(@[[ join(",", (":") x $dim) ]]), intent(in) :: a
    write(*,*) "A @[[$dim]] dimensional array of shape", shape(a)
  end subroutine

!PERL-END }

end module perl_test_mod

!make NO_CONFIG=1
program perl_test
  use perl_test_mod
  use abort
  implicit none

  real :: a(4), b(4,4), c(4,4,4)

  call array_info(a)
  call array_info(b)
  call array_info(c)

!PERL-START my @number = (1, 2, 3);
!PERL       my @letter = ("a", "b", "c");
!PERL       for (my $i = 0; $i < scalar(@number); $i++) {
!PERL           my $n = $number[$i];
!PERL           my $l = $letter[$i];
  write(*,*) @[[$n]], "@[[$l]]"
!PERL-END }

end program
