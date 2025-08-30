##############################################################
# dependencylist.pl generated file, Mon Aug 25 15:53:59 2025 #
##############################################################

programs = copy dump_test env_test eos_python perl_test test_diff test_hdf5_io test_interpolation test_newton_raphson test_stringutils test_timings vertex

# program copy in copy.F90
mainfile_copy = copy.F90

copy:
	@$(MAKE) --no-print-directory MAIN=copy CONFIG=make.inc.vertex

copy-info:
	@$(MAKE) --no-print-directory MAIN=copy info CONFIG=make.inc.vertex

# program dump_test in dump_test.F90
mainfile_dump_test = dump_test.F90

dump_test:
	@$(MAKE) --no-print-directory MAIN=dump_test 

dump_test-info:
	@$(MAKE) --no-print-directory MAIN=dump_test info 

# program env_test in env_test.F90
mainfile_env_test = env_test.F90

env_test:
	@$(MAKE) --no-print-directory MAIN=env_test NO_CONFIG=1 DEBUG=HIGH

env_test-info:
	@$(MAKE) --no-print-directory MAIN=env_test info NO_CONFIG=1 DEBUG=HIGH

# program eos_python in eos_wrap.F90
mainfile_eos_python = eos_wrap.F90

eos_python:
	@$(MAKE) --no-print-directory MAIN=eos_python 

eos_python-info:
	@$(MAKE) --no-print-directory MAIN=eos_python info 

# program perl_test in perl_test.F90
mainfile_perl_test = perl_test.F90

perl_test:
	@$(MAKE) --no-print-directory MAIN=perl_test NO_CONFIG=1

perl_test-info:
	@$(MAKE) --no-print-directory MAIN=perl_test info NO_CONFIG=1

# program test_diff in derivative.F90
mainfile_test_diff = derivative.F90

test_diff:
	@$(MAKE) --no-print-directory MAIN=test_diff NO_CONFIG=1

test_diff-info:
	@$(MAKE) --no-print-directory MAIN=test_diff info NO_CONFIG=1

# program test_hdf5_io in mod_dataformat_hdf5_test.F90
mainfile_test_hdf5_io = mod_dataformat_hdf5_test.F90

test_hdf5_io:
	@$(MAKE) --no-print-directory MAIN=test_hdf5_io NO_CONFIG=1 EXTRA_CPPFLAGS="-DUNIT_TESTS -DMPI_HYDRO" DEBUG=HIGH ENABLE=MPI

test_hdf5_io-info:
	@$(MAKE) --no-print-directory MAIN=test_hdf5_io info NO_CONFIG=1 EXTRA_CPPFLAGS="-DUNIT_TESTS -DMPI_HYDRO" DEBUG=HIGH ENABLE=MPI

# program test_interpolation in interpolation.F90
mainfile_test_interpolation = interpolation.F90

test_interpolation:
	@$(MAKE) --no-print-directory MAIN=test_interpolation NO_CONFIG=1 EXTRA_CPPFLAGS=-DUNIT_TESTS

test_interpolation-info:
	@$(MAKE) --no-print-directory MAIN=test_interpolation info NO_CONFIG=1 EXTRA_CPPFLAGS=-DUNIT_TESTS

# program test_newton_raphson in newton_raphson.F90
mainfile_test_newton_raphson = newton_raphson.F90

test_newton_raphson:
	@$(MAKE) --no-print-directory MAIN=test_newton_raphson NO_CONFIG=1 EXTRA_CPPFLAGS=-DUNIT_TESTS

test_newton_raphson-info:
	@$(MAKE) --no-print-directory MAIN=test_newton_raphson info NO_CONFIG=1 EXTRA_CPPFLAGS=-DUNIT_TESTS

# program test_stringutils in strings.F90
mainfile_test_stringutils = strings.F90

test_stringutils:
	@$(MAKE) --no-print-directory MAIN=test_stringutils NO_CONFIG=1

test_stringutils-info:
	@$(MAKE) --no-print-directory MAIN=test_stringutils info NO_CONFIG=1

# program test_timings in new_timings.F90
mainfile_test_timings = new_timings.F90

test_timings:
	@$(MAKE) --no-print-directory MAIN=test_timings NO_CONFIG=1 DEBUG=HIGH

test_timings-info:
	@$(MAKE) --no-print-directory MAIN=test_timings info NO_CONFIG=1 DEBUG=HIGH

# program vertex in vertex.F90
mainfile_vertex = vertex.F90

vertex:
	@$(MAKE) --no-print-directory MAIN=vertex 

vertex-info:
	@$(MAKE) --no-print-directory MAIN=vertex info 


$(ppdir)/burn.f90: \
	src/burning/flashing.X90

$(ppdir)/high_density_eos.f90: \
	src/eos/high_density_eos.trilinear.X90

# vi: syntax=make
