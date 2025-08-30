## *** MAKEFILE TARGETS ***
##
##  - no target - or "all"
##	Compile "vertex"
##
##  "help"
##	Gives this help
##
##  "help-features"
##	Some information about all this compilation feature stuff
##
##  "info"
##	List the compiler and build flags as well as the preprocessor
##	directives that would be used
##
##  "clean"
##	rm -rf build/
##
##  Any name of a Fortran "program" statement
##	Collects all required fortran source files for this starting
##	point, preprocesses and compiles them into an executable
##
##	Actually just a synonym for
##
##		make MAIN=programname
##
##  "list-programs"
##	Print out a list of programs in the source files
##
##  "preprocess"
##	Preprocess all fortran files
##
##  "$(ppdir)/somefile.f90" or "$(ppdir)/somefile.f"
##	Run the preprocessor with the current configuration to produce
##	only a specific file (due the the way dependencies are tracked,
##	all files will actually be preprocessed and you will get the
##	message, that the target is already up to date)
##
##  "$(builddir)/somefile.o"
##	Build only a specifiy object file with the current configuration
##
##  "defined-flags"
##	List all preprocesser flags that are "defined", i.e. specified with -D
##
##  "print_executable_name"
##	Echo the filename that make all would build, for scripting
##
##  "check-preprocessor"
##	Run a (crude!) perl script to enforce a minimally consistent use of
##	preprocessor statements. Part of the compiletest
##
##  "check-cycle"
##	Run a script to check, wheter the dependency graph generated via fortran
##	"use" statements is acyclic, which is necessary for compilation. Part of
##	the compiletest
##
##  "check-no-tabs"
##	Run a script to ensure that there are no TAB characters present in the
##	fortran source file. Part of the compiletest
##
##  "test"
##	Run all check-* targets on the default configuration and try to compile
##	every configuration found in configurations/.
##
##  A filename in configurations
##	Run a selective compile test on a make.inc.vertex like file in configurations/
##
##  "run-on-pp"
##	Preprocess all files with the current configuration and run a tool
##	specified with the P parameter on those files, to be used like
##
##		make run-on-pp P="program {}"
##
##      The "{}" is replaced with the preprocessed fortran files.
##
##  "doc"
##	Use doxygen to generate documentaion in doc/
##
##  "tags"
##	Generate a tags file for use with your editor
##
##  "forcheck"
##	Run forcheck on preprocessed source code
##
##  "why-need"
##      Query why a certain file is necessary, use it as
##
##              make why-need FILE=eos.f90
##
##      to get the shortest "use" path from eos.f90 to the main fortran file
##      with the "program" statement
##
## *** GENERAL REMARKS ***
##
## Put Fortran or C-files within a directory in VPATH, the dependencies between
## Fortran files via "use module" statements and for both Fortran and C-Files
## via "#include" directives are honored and extracted with the script
## dependencylist.pl in tools/
##
##
## *** MAKEFILE OPTIONS ***
##
## NOSILENT=1
##   Set to 1 to echo compilation commands to the screen
##
## MAIN = string
##   Set to the name of the Fortran "program" that you want as main. All
##   necessary files for compilation are gathered via the "use" statements that
##   start in the file containing this "program".
##
## ENABLE = list
##   Set to a list of compilation "features" (see below) that you want to add to
##   the default ones
##
## DISABLE = list
##   Specify a list of default features not to use for this compilation
##
## WITH_FEATURE = {0, 1}
##   Optional syntax for disabling/enabling individual features
##
## DEBUG = {LOW, MEDIUM, HIGH}
##   - DEBUG=LOW:
##	enable the feature DEBUG_SYMBOLS
##
##   - DEBUG=MEDIUM
##	enable DEBUG_SYMBOLS and NO_OPT
##
##   - DEBUG=HIGH
##	enable DEBUG_SYMBOLS, NO_OPT, and DEBUG_CHECKS
##
## EXTRA_CPPFLAGS
##   Specify a list of flags for the preprocessor that get appended to CPPFLAGS
##
##
#F *** WHAT IS A COMPILATION "FEATURE"? ***
#F
#F The architecture makefiles can specify different modes of compilation in a
#F modular fashion by defining certain "features". A feature is just a set of
#F specially named variables which act on the compilation commands if the
#F feature is active.
#F
#F For a feature to work its name must be registered in the main Makefile in
#F the known_features variable.
#F
known_features = HDF5 \
		 HDF5_SHARED \
		 OPENMP \
		 OPT \
		 NO_OPT \
		 DEBUG_CHECKS \
		 DEBUG_SYMBOLS \
		 MPI \
		 VENDOR_VECLIBS \
		 VENDOR_LAPACK \
		 GSL \
		 STRICT_FP \
		 PIC \
		 MARMOT \
		 THREAD_CHECKER \
		 PERFLIB \
		 SCALASCA \
		 TRACE_ANALYSER \
		 MPITRACE \
		 CALL_TRACE \
		 EXPORT_DYNAMIC	\
		 CUDAF \
		 CUDAC
export known_features
#F
#F *** A SHORT DESCRIPTION OF THE DEFAULT FEATURES ***
#F
#F  HDF5                link with HDF5 library, active by default
#F  OPENMP              compile with OpenMP support, active by default
#F  OPT                 optimize code, active by default
#F
#F  NO_OPT              do not optimize code
#F  DEBUG_CHECKS        compile with aggressive debug flags, such as array bounds
#F                        checking and floating-point traps
#F  DEBUG_SYMBOLS       compile with debug symbols
#F  MPI                 link with MPI libraries
#F  VENDOR_VECLIBS      use fast vendor libraries for math functions
#F  VENDOR_LAPACK       use fast vendor linear algebra routines
#F  STRICT_FP           If possible, enforce strict IEEE conforming floating
#F                        point operations
#F  PIC                 Generate position independent code (for use in a library)
#F  MARMOT              Use marmot compiler frontend
#F  PERFLIB             Link with perflib
#F  SCALASCA            Use Scalasca
#F  MPITRACE            Link with mpitrace lib
#F  CALL_TRACE          Print a debug message on every function call and exit
#F  EXPORT_DYNAMIC      Instruct linker to add all symbols to the dynamic symbol table
#F                        This is necessary for any run-time symbol resolution such as
#F                        used in CALL_TRACE or via -DBACKTRACE_ON_ABORT
#F
#F
#F Architecture sub makefiles may define their own, private features by adding
#F to known_features within the sub makefile itself.
#F
#F
#F *** HOW FEATURES WORK INTERNALLY ***
#F
#F For a feature called "FEATURENAME" you can set any of these variables
#F altering the behaviour of the compilation:
#F
#F   FEATURENAME_CPPFLAGS
#F   FEATURENAME_FFLAGS
#F   FEATURENAME_CFLAGS
#F   FEATURENAME_LDFLAGS
#F
#F   FEATURENAME_FC
#F   FEATURENAME_CC
#F   FEATURENAME_LD
#F
#F   FEATURENAME_FC_WRAPPER
#F   FEATURENAME_CC_WRAPPER
#F   FEATURENAME_LD_WRAPPER
#F
#F   FEATURENAME_HDF5_HOME
#F
#F If FEATURENAME is active the appropriate values get appended to the CCPFLAGS,
#F FFLAGS, CFLAGS, and LDFLAGS variables. If more than one feature adds values,
#F they are appended in lexicographical order of the feature's names.
#F
#F The default values of FC, CC, LD, and HDF5_HOME are replaced(!) should a
#F feature specify those, therefore no more than one active feature may specify
#F a compiler, linker, or HDF5_HOME; an error is raised in case there is more
#F than one definiton active.
#F
#F The _WRAPPER variables are prepended (again in lexicographical order of the
#F feature's names) to the current value of FC, CC, or LD to allow specifying
#F some wrapper programs (e.g. scalasca).
#F
#F In order to declare possible conflicts, one may specify
#F
#F   FEATURENAME_CONFLICTSWITH
#F
#F and set it to a list of other features which should not be used together with
#F feature FEATURENAME, resulting in an error when trying to invoke make. Also,
#F one may specify
#F
#F   FEATURENAME_REQUIRES
#F
#F and set it to a list of other features which have to be specified along with
#F this feature, an error is raised otherwise


# Actual Makefile follows

#SHELL=#/bin/bash

# No built-in suffix rules
.SUFFIXES:

# macro for verbose or short echoing
ifeq ($(NOSILENT),1)
  echo=echo "$(strip $1)"; $1
else
  echo=$(if $2,echo "$2";,) $1 || (echo -en "\n\nError running\n  $1\n\n" 1>&2; exit 1;)
endif

CPP = cpp -traditional -P -Wall -Werror

example_configs = $(notdir $(wildcard configurations/*))
phonies += $(programs)

VPATH = src/burning \
        src/burning/modules \
        src/c-code \
        src/c-code/debug \
        src/c-code/setup \
        src/c-code/memory \
        src/coconut \
        src/coconut/modules \
        src/debug \
        src/eos \
        src/eos/data \
        src/hydro \
        src/hydro/modules \
        src/io \
        src/io/modules \
        src/libs \
        src/parallelisation \
        src/parallelisation/modules \
        src/middleware \
        src/startup \
        src/system \
        src/system/modules \
        src/transport \
        src/transport/modules \
        src/util \
        src/python \
        src/tests

ifndef MAIN
  MAIN=vertex
endif

id=$(PREFIX)$(MAIN)-$(ARCH)
sourcedep = source.dep.mk
ppdep = build/$(id)/pp.dep.mk
depgraph = build/$(id)/deps.dot

# wheter to load the dependencies and config files
no_config_targets = clean \
		    help \
		    help-features \
		    doc-clean \
		    test \
		    check-preprocessor \
		    check-no-tabs \
		    %.so \
		    %.so.info \
		    $(example_configs) \
		    list-cherries

no_dep_targets = clean \
		 help \
		 help-features \
		 doc-clean \
		 test \
		 check-preprocessor \
		 check-no-tabs \
		 info \
		 $(programs) \
		 $(addsuffix -info,$(programs)) \
		 list-programs \
		 defined-flags \
		 %.so \
		 %.so.info \
		 $(example_configs) \
		 list-cherries

no_arch_targets = clean \
		  help \
		  help-features \
		  doc-clean \
		  check-no-tabs \
		  list-cherries

NO_CONFIG=$(if $(filter $(no_config_targets),$(MAKECMDGOALS)),1,0)
NO_DEP=$(if $(filter $(no_dep_targets),$(MAKECMDGOALS)),1,0)
NO_ARCH=$(if $(filter $(no_arch_targets),$(MAKECMDGOALS)),1,0)

# Vertex executable is named specially for historic reasons
exec=$(if $(filter vertex,$(MAIN)),build/$(id)/$(EXE_DEFAULT).x.$(HOST),build/$(id)/$(MAIN).x.$(HOST))

# the default since first target
default_target: all

include $(sourcedep)

# Check wheter $(MAIN) is a valid program
define newline


endef
ifdef programs
nothing := $(if $(filter $(MAIN),$(programs)),,$(error No such program "$(MAIN)".$(newline)Valid programs are:$(newline)$(addprefix $(newline),$(sort $(programs)))$(newline)$(newline)Make sure there is a "program $(MAIN)" statement in the sources))
endif

ifneq ($(NO_CONFIG),1)
# load make.inc.$(MAIN) file when $(sourcedep) was made and loaded
# i.e. "programs" is defined
  ifdef programs
    ifndef CONFIG
      CONFIG=make.inc.$(MAIN)
    endif
    include $(CONFIG)
    CPPFLAGS += -DHAVE_CONFIG
  endif
endif

ppdir=build/$(id)/pp
builddir=build/$(id)/objects

# Create the pp/, build/ directories in case they do not exist yet
ifneq ($(NO_ARCH),1)
  nothing := $(shell mkdir -p $(ppdir) $(builddir))
endif

# List of files
sources=$(notdir $(shell find $(VPATH) -maxdepth 1 -name "*.F90" -o -name "*.F" -o -name "*.c" -o -name "*.cu" -o -name "*.CUF"))
f90_ppsources = $(patsubst %.F90, $(ppdir)/%.f90, $(filter %.F90, $(sources)))
f_ppsources   = $(patsubst %.F,   $(ppdir)/%.f,   $(filter %.F,   $(sources)))
cudaf_ppsources   = $(patsubst %.CUF,   $(ppdir)/%.cuf,   $(filter %.CUF,   $(sources)))
ppsources = $(f90_ppsources) $(f_ppsources) $(cudaf_ppsources)
c_objects    = $(addprefix $(builddir)/, $(patsubst %.c, %.o, $(filter %.c, $(sources))))
cudaf_objects = $(if $(WITH_CUDAF).$(addprefix $(builddir)/, $(patsubst %.cuf, %.o, $(filter %.cuf, $(sources)))),)
cudac_objects = $(if $(WITH_CUDAC).$(addprefix $(builddir)/, $(patsubst %.cu, %.o, $(filter %.cu, $(sources)))),)

text_objects = $(addprefix $(builddir)/, version_diff.o version_log.o version_commit_id.o)
ifdef CONFIG
  text_objects += $(builddir)/$(notdir $(CONFIG)).o
endif

ifneq ($(NO_ARCH),1)
  ifndef ARCH
    $(error You have to specify a file within Makefile.arch/ via the ARCH Parameter or Environment variable)
  endif

  # macro for use within Makefile.arch/$(ARCH)
  #
  # Usage:   $(call unless_feature_set,FEATURENAMES,STRING)
  #
  # emits STRING unless one of the features in FEATURENAMES is set
  # e.g. to omit
  #
  dummy_tag = __VERTEX_DUMMY__
  unless_feature_set = $(if $(strip $(filter $1,$(FEATURES))),$(dummy_tag),$2)

  # Get options from sub makefiles
  include Makefile.arch/INCLUDE
  include Makefile.arch/$(ARCH)

  # Features that are on by default:
  # - openmp
  # - optimization
  # - hdf5
  DEFAULT_FEATURES = OPENMP OPT HDF5

  FEATURES += $(DEFAULT_FEATURES)

  enable = $(ENABLE)
  disable = $(DISABLE)

  # Enable/disable features by WITH_FEATURENAME=1 or 0
  # sanity check
  nothing := $(foreach feature,$(known_features),$(if $(filter-out 1 0,$(WITH_$(feature))),$(error WITH_$(feature) must be 0 or 1),))
  # add them to the ENABLE, DISABLE vars
  enable  += $(foreach feature,$(known_features),$(if $(findstring 1,$(WITH_$(feature))),$(feature),))
  disable += $(foreach feature,$(known_features),$(if $(findstring 0,$(WITH_$(feature))),$(feature),))

  ifdef DEBUG
    ifeq ($(DEBUG),HIGH)
      FEATURES += DEBUG_SYMBOLS DEBUG_CHECKS
      disable += OPT
    endif
    ifeq ($(DEBUG),MEDIUM)
      FEATURES += DEBUG_SYMBOLS
      disable += OPT
    endif
    ifeq ($(DEBUG),LOW)
      FEATURES += DEBUG_SYMBOLS
    endif
    ifneq (,$(filter-out LOW MEDIUM HIGH,$(DEBUG)))
      $(error Unknown entry in DEBUG: $(DEBUG))
    endif
  endif

  ifneq (,$(filter-out $(known_features),$(enable)))
    $(error Unknown entries in ENABLE: $(filter-out $(known_features),$(enable)))
  endif
  FEATURES += $(enable)

  ifneq (,$(filter-out $(known_features),$(disable)))
    $(error Unknown entries in DISABLE: $(filter-out $(known_features),$(disable)))
  endif

  # You can disable features by calling make with
  # a DISABLE variable, e.g.
  #
  # make DISABLE="MPI OPT" target
  #
  # filter-out disabled features
  FEATURES := $(filter-out $(disable),$(FEATURES))

  # if no OPT use NO_OPT
  ifeq (,$(filter OPT,$(FEATURES)))
    FEATURES += NO_OPT
  endif

  # remove duplicates and sort
  FEATURES := $(sort $(FEATURES))
  ifneq (,$(filter-out $(known_features),$(FEATURES)))
    $(error Unknown entries in FEATURES: $(filter-out $(known_features),$(FEATURES)))
  endif

  # Test if some feature was wanted but does not specify anything
  # in this architecture makefile
  possible_feature_flags = CPPFLAGS FFLAGS CFLAGS CUDACFLAGS CUDAFFLAGS LDFLAGS FC CC CUDAFC CUDACC LD HDF5_HOME FC_WRAPPER CC_WRAPPER LD_WRAPPER
  nothing := $(foreach feature,$(FEATURES),$(if $(strip $(foreach flag,$(possible_feature_flags),$($(feature)_$(flag)))),,$(error Feature $(feature) is undefined for ARCH=$(ARCH))))

  # Test if a feature conflicts with another
  nothing :== $(foreach feature,$(FEATURES),$(if $(filter $($(feature)_CONFLICTSWITH),$(FEATURES)),$(error Feature $(feature) conflicts with other enabled feature(s): $(filter $($(feature)_CONFLICTSWITH),$(FEATURES)))))

  # Test if all requirements are met
  nothing :== $(foreach feature,$(FEATURES),$(if $(strip $(filter-out $(FEATURES),$($(feature)_REQUIRES))),$(error Feature $(feature) needs the feature(s): $(filter-out $(FEATURES),$($(feature)_REQUIRES)))))

  # Special variable for MPI Compilation
  ifneq (,$(findstring MPI,$(FEATURES)))
    ifeq ($(USE_MPIF_H),1)
      CPPFLAGS += -DUSE_MPIF_H
    endif
  endif

  # Append feature flags to actual flags
  feature_flags = $(foreach feature,$(FEATURES),$(filter-out $(dummy_tag),$($(feature)_$(1))))

  CPPFLAGS += $(call feature_flags,CPPFLAGS)
  CPPFLAGS += $(foreach feature,$(FEATURES),-DWITH_$(feature))

  # Set WITH_FEATURENAME=1 for all enabled features, clear for all unset
  $(foreach feature,$(known_features),$(eval $(if $(filter $(feature),$(FEATURES)),WITH_$(feature):=1,WITH_$(feature):=)))

  FFLAGS   += $(call feature_flags,FFLAGS)
  CFLAGS   += $(call feature_flags,CFLAGS)
  CUDACFLAGS+= $(call feature_flags,CUDACFLAGS)
  CUDAFFLAGS+= $(call feature_flags,CUDAFFLAGS)
  LDFLAGS  += $(call feature_flags,LDFLAGS)

  # See if we have to replace the default FC, CC, LD, and HDF5_HOME. Emit an
  # error if there is more than one non-default value for those defined
  comma_separated_list = $(shell foo=($1); IFS=","; echo "$${foo[*]}";)
  select_choice = $(strip $(call select_choice_,$1,$(call feature_flags,$1)))
  select_choice_ = $(if $(filter-out 0 1,$(words $2)),$(error More than one $1 specified: $(call comma_separated_list,$2)), \
                     $(if $(filter 0,$(words $2)),$($1),$2))
  FC := $(call feature_flags,FC_WRAPPER) $(call select_choice,FC)
  CC := $(call feature_flags,CC_WRAPPER) $(call select_choice,CC)
  CUDAFC := $(call feature_flags,CUDAFC_WRAPPER) $(call select_choice,CUDAFC)
  CUDACC := $(call feature_flags,CUDACC_WRAPPER) $(call select_choice,CUDACC)
  LD := $(call feature_flags,LD_WRAPPER) $(call select_choice,LD)
  HDF5_HOME := $(call select_choice,HDF5_HOME)

  # Set a default routine for the fortran module file name translation
  FORTRAN_MOD_FILE_NAME := $(if $(FORTRAN_MOD_FILE_NAME),$(FORTRAN_MOD_FILE_NAME),'sub { my $$mod = shift; return lc($$mod) . ".mod"; }')
endif

ifneq ($(NO_DEP),1)
# Generate Fortran module dependencies when $(sourcedep) was made and loaded
  ifdef programs
    include $(ppdep)
  endif

  # Location of tag files for the various flag variables
  flag_dep_prefix = $(ppdir)/__FLAG_FILE_
  flag_dep = $(addprefix $(flag_dep_prefix),$(1))

  # this needs secondary expansion
  flag_dep_add_remove_replace = $(foreach flag,$(2),$(foreach mod,ADD REMOVE REPLACE,$$(if $$($(flag)_$(mod)_$$*$(strip $(1))),$$(call flag_dep, $(flag)_$(mod)_$$*$(strip $(1))),)))

  # Update dependency files for changed flags
  # The file for FLAG is only modified if the content changed, thus
  # triggering a recompilation all rules with the prerequisite $(call flag_dep, FLAG)
  ifeq ($(DEBUG_FLAG_DEPS),1)
    update_flag_deps = $(shell echo "$($(1))" > "$(call flag_dep,$(1))__"; if [ ! -e "$(call flag_dep,$(1))" ] || ! cmp -s "$(call flag_dep,$(1))__" "$(call flag_dep,$(1))"; then echo -e "Updating flag dependency for $(1) in $(id)\n$$(diff -u $(call flag_dep,$(1))__ $(call flag_dep,$(1)))\n" 1>&2; mv "$(call flag_dep,$(1))__" "$(call flag_dep,$(1))"; else rm "$(call flag_dep,$(1))__"; echo "Not updating flag dependency for $(1) in $(id)" 1>&2; fi; )
  else
    update_flag_deps = $(shell echo "$($(1))" > "$(call flag_dep,$(1))__"; if [ ! -e "$(call flag_dep,$(1))" ] || ! cmp -s "$(call flag_dep,$(1))__" "$(call flag_dep,$(1))"; then mv "$(call flag_dep,$(1))__" "$(call flag_dep,$(1))"; else rm "$(call flag_dep,$(1))__"; fi; )
  endif

  flags = $(possible_feature_flags) CPP FORTRAN_MOD_FILE_NAME FFLAGS_FREE FFLAGS_FIXED

  # foreach defined generic flag, create a flag-file
  nothing := $(foreach flag,$(flags),$(if $(filter-out undefined,$(origin $(flag))),$(call update_flag_deps,$(flag)),))

  # foreach defined file specific flag, create a flag file
  nothing := $(foreach file,$(sources) $(patsubst %.F90, %.f90, $(filter %.F90, $(sources))),$(foreach flag,$(flags),$(foreach mod,ADD REMOVE REPLACE,$(if $($(flag)_$(mod)_$(file)),$(call update_flag_deps,$(flag)_$(mod)_$(file))))))

  # here it gets complicated:
  # In case a file-specific flag has been previously defined but is now unset,
  # we also have to trigger a recompilation of the affected file.
  # This is accomplished by "touch"ing the file itself.
  touch_source_file = $(shell find $(VPATH) $(ppdir) $(builddir) -maxdepth 1 -name "$(1)" -exec touch {} \;)
  filename_of_flag_variable = $(strip $(foreach flag,$(flags),$(foreach mod,ADD REMOVE REPLACE,$(if $(findstring $(flag)_$(mod),$(1)),$(patsubst $(flag)_$(mod)_%,%,$(1)),))))
  nothing := $(foreach file,$(wildcard $(call flag_dep,*)),$(if $(filter-out undefined,$(origin $(patsubst $(flag_dep_prefix)%,%,$(file)))),,$(call touch_source_file,$(call filename_of_flag_variable, $(patsubst $(flag_dep_prefix)%,%,$(file))))$(shell rm $(file))))
endif

mainfile = $(mainfile_$(MAIN))

# global preprocessor flags
CPPFLAGS += -DPROGRAM_$(MAIN) $(EXTRA_CPPFLAGS) -UDOUBLE_PRECISION_EOS
global_include=src/macros/header.h

# perl postprocessor
post_cpp_f90=src/macros/post_cpp.f90.pl
post_cpp_f77=src/macros/post_cpp.f77.pl

# target to preprocess all source files
preprocess: $(ppsources)

configfiles=$(CONFIG) Makefile.arch/$(ARCH) Makefile

# Executable file
all: $(exec)

$(exec): $(fortran_objects) $(c_objects) $(cudac_objects) $(cudaf_objects) $(text_objects) $(builddir)/compile_time.o $(builddir)/version.o $(call flag_dep,LDFLAGS)
	@$(call echo,$(LD) $(filter-out $(call flag_dep,LDFLAGS),$^) $(LDFLAGS) -o $@,LD  $(<F))
	@echo
	@echo "Compiled $(exec)"

print_executable_name:
	@echo $(exec)

# Dependencies
$(sourcedep): $(sources) tools/dependencies_unprocessed.pl tools/dependency.pm $(call flag_dep,CPPFLAGS)
	@SOURCES="$(filter src/%,$^)" \
	 INCLUDEDIRS=$(subst -I,,$(filter -I%,$(CPPFLAGS))) \
	 tools/dependencies_unprocessed.pl > $@ || (rm $@ && false)

$(ppdep): $(ppsources) $(sourcedep) tools/dependencies_preprocessed.pl tools/dependency.pm $(call flag_dep,FORTRAN_MOD_FILE_NAME)
	@tools/dependencies_preprocessed.pl $(mainfile) $(ppdir) FORTRAN_MOD_FILE_NAME=$(FORTRAN_MOD_FILE_NAME) > $@ || (rm $@ && false)

$(depgraph): $(ppsources) $(sources) tools/dependency.pm tools/dependencygraph.pl $(call flage_dep,CPPFLAGS FORTRAN_MOD_FILE_NAME)
	@SOURCES="$(filter src/%,$^)" \
	 PPDIR="$(ppdir)" \
	 INCLUDEDIRS="$(subst -I,,$(filter -I%,$(CPPFLAGS)))" \
	 FORTRAN_MOD_FILE_NAME=$(FORTRAN_MOD_FILE_NAME) \
	 tools/dependencygraph.pl $(MAIN) > $@ || (rm $@ && false)

# List of programs
list-programs: $(sourcedep)
	@echo Valid programs are:
	@echo
	@$(foreach prog,$(sort $(programs)),echo -en "- \"$(prog)\" in $(mainfile_$(prog))\n";)
	@echo

# general rule set:
#
# Preprocessing
#
define preprocess
  if [ -e "$@" ] ; then chmod +w $@; fi
  $(call echo,$(CPP) $(call cppflags,$(<F)) $< -o $@,CPP $(<F))
endef

define preprocess_imacros
  if [ -e "$@" ] ; then chmod +w $@; fi
  $(call echo,$(CPP) $(call cppflags,$(<F)) -imacros $(global_include) $< -o $@,CPP $(<F))
endef

.SECONDEXPANSION:
$(ppdir)/error.f90 : error.F90 $(call flag_dep, CPP CPPFLAGS) $(call flag_dep_add_remove_replace, .F90, CPPFLAGS)
	@$(preprocess)
	@chmod -w $@

$(ppdir)/abort.f90 : abort.F90 $(call flag_dep, CPP CPPFLAGS) $(call flag_dep_add_remove_replace, .F90, CPPFLAGS)
	@$(preprocess)
	@chmod -w $@

$(ppdir)/%.f90 : %.F90 $(global_include) $(call flag_dep, CPP CPPFLAGS) $(call flag_dep_add_remove_replace, .F90, CPPFLAGS)
	@$(preprocess_imacros)
	@chmod -w $@

$(ppdir)/%.cuf : %.CUF $(global_include) $(call flag_dep, CPP CPPFLAGS) $(call flag_dep_add_remove_replace, .CUF, CPPFLAGS)
	@$(preprocess_imacros)
	@chmod -w $@

$(ppdir)/%.f : %.F $(global_include) $(call flag_dep, CPP CPPFLAGS) $(call flag_dep_add_remove_replace, .F, CPPFLAGS)
	@$(preprocess_imacros)
	@chmod -w $@

# Implement {C,F,CUDAC,CUDAF}FLAGS_{ADD,REMOVE,REPLACE}_somefile logic in one place
define flags
  $(if $($(1)FLAGS_REPLACE_$(2)),$(if $($(1)FLAGS_ADD_$(2))$($(1)FLAGS_REMOVE_$(2)),$(error You cannot specify $(1)FLAGS_REPLACE_$(2) together with any of $(1)FLAGS_ADD_$(2) or $(1)FLAGS_REMOVE_$(2)),)$($(1)FLAGS_REPLACE_$(2)),$(filter-out $($(1)FLAGS_REMOVE_$(2)),$($(1)FLAGS)) $($(1)FLAGS_ADD_$(2)))
endef

define cppflags
  $(call flags,CPP,$(1))
endef

define fflags
  $(call flags,F,$(1))
endef

define cflags
  $(call flags,C,$(1))
endef

define cudafflags
  $(call flags,CUDAFFLAGS,$(1))
endef

define cudacflags
  $(call flags,CUDACFLAGS,$(1))
endef

# Compilation
#
#  Fortran files
$(builddir)/%.o : $(builddir)/%.f90 $(call flag_dep, FC FFLAGS_FREE FFLAGS) $(call flag_dep_add_remove_replace, .f90, FFLAGS)
	@cd $(builddir) && $(call echo,$(FC) $(call fflags,$(<F)) $(FFLAGS_FREE) -c $(<F) -o $*.o,FC  $(<F))
	@$(if $(file_$(<F)_modules),touch -c -r $(builddir)/$*.o $(addprefix $(builddir)/,$(file_$(<F)_modules)),)

$(builddir)/%.o : $(builddir)/%.f $(call flag_dep, FC FFLAGS_FIXED FFLAGS) $(call flag_dep_add_remove_replace, .f, FFLAGS)
	@cd $(builddir) && $(call echo,$(FC) $(call fflags,$(<F)) $(FFLAGS_FIXED) -c $(<F) -o $*.o,FC  $(<F))
	@$(if $(file_$(<F)_modules),touch -c -r $(builddir)/$*.o $(addprefix $(builddir)/,$(file_$(<F)_modules)),)

# It is nice to have a builddir whith object and
# source files together to make some debuggers happy
_space_mod = $(if $(filter-out 0,$(NO_ARCH)),"You_have_to_load_an_architecture_makefile_to_build_modules",$(shell perl -e 'print (('$(FORTRAN_MOD_FILE_NAME)')->(" "));'))
_mod_post = $(lastword $(_space_mod))
_mod_pre = $(patsubst %$(_mod_post),%,$(_space_mod))
.PRECIOUS: $(builddir)/%.f90 $(builddir)/%.f $(builddir)/%.c $(ppdir)/%.f90 $(ppdir)/%.f $(ppdir)/%.c $(builddir)/$(_mod_pre)%$(_mod_post)

# to make SCALASCA happy, insert a comment as first line
$(builddir)/%.f : $(ppdir)/%.f $(post_cpp_f77)
	@echo '! $(@F)' > $@
	@($(post_cpp_f77) $< >> $@) || (rm $@ && false)

# to make SCALASCA happy, insert a comment as first line
$(builddir)/%.f90 : $(ppdir)/%.f90 $(post_cpp_f90)
	@echo '! $(@F)' > $@
	@($(post_cpp_f90) $< >> $@) || (rm $@ && false)

$(builddir)/%.cuf : $(ppdir)/%.cuf $(post_cpp_f90)
	@echo '! $(@F)' > $@ # to make SCALASCA happy, insert a comment as first line
	@($(post_cpp_f90) $< >> $@) || (rm $@ && false)

$(builddir)/%.c : $(ppdir)/%.c
	@cd $(builddir) && ln -sf ../pp/$(<F) .

# cuf-files from pp
$(builddir)/%.o : $(builddir)/%.cuf $(call flag_dep, CUDAFC CUDAFFLAGS) $(call flag_dep_add_remove_replace, .cuf, CUDAFFLAGS)
	@cd $(builddir) && $(call echo,$(CUDAFC) $(call cudafflags,$(<F)) -c $(<F) -o $*.o,CUF $(<F))
	@$(if $(file_$(<F)_modules),touch -c -r $(builddir)/$*.o $(addprefix $(builddir)/,$(file_$(<F)_modules)),)

#  cu-files from source
$(builddir)/%.o : %.cu $(call flag_dep, CUDACC CPPFLAGS CUDACFLAGS) $(call flag_dep_add_remove_replace, .cu, CPPFLAGS CUDACFLAGS)
	@$(call echo,$(CUDACC) $(call cppflags,$(<F)) $(call cudacflags,$(<F)) $< -c -o $@,CUC $(<F))

#  c-files from source
$(builddir)/%.o : %.c $(call flag_dep, CC CPPFLAGS CFLAGS) $(call flag_dep_add_remove_replace, .c, CPPFLAGS CFLAGS)
	@$(call echo,$(CC) $(call cppflags,$(<F)) $(call cflags,$(<F)) $< -c -o $@,CC  $(<F))

# sometimes one wants to look at a preprocessed c-file
$(builddir)/%.i : %.c $(call flag_dep, CC CPPFLAGS CFLAGS) $(call flag_dep_add_remove_replace, .c, CPPFLAGS CFLAGS)
	@$(call echo,$(CC) $(call cppflags,$(<F)) $(call cflags,$(<F)) $< -E > $@,CC  $(<F))

# c-files from pp
$(builddir)/%.o : $(ppdir)/%.c $(call flag_dep, CC CPPFLAGS CFLAGS) $(call flag_dep_add_remove_replace, .c, CPPFLAGS CFLAGS)
	@$(call echo,$(CC) $(call cppflags,$(<F)) $(call cflags,$(<F)) $< -c -o $@,CC  $(<F))

# Include a file's content into an object in a portable way
define bin_to_c_file
  (echo "#include <stdlib.h>"; \
   echo "const char $(1)[] = {"; \
   od -vtxC $(2) | sed -e "s/^[0-9]*//" -e s"/ \([0-9a-f][0-9a-f]\)/0x\1,/g" -e"\$$d" ; \
   echo "0x00};") > $(3); \
  echo "const size_t $(1)_length = `wc -c < $(2)`;" >> $(3)
endef

# Generated files with info texts
$(ppdir)/version_diff.c : $(sources) $(configfiles) 
#	@if hash git 2>/dev/null; then \
#		echo "git diff blessed" > $(ppdir)/version_diff; \
#		git diff blessed >> $(ppdir)/version_diff; \
#	else \

		sleep 10; \
		echo 1>&2; \
		echo "WARNING! No git installed on build host. Not saving any diffs!" 1>&2; \
		echo 1>&2; \
		echo "No git installed on build host" > $(ppdir)/version_diff; \
#	fi;
	@$(call bin_to_c_file, _version_diff, $(ppdir)/version_diff, $(ppdir)/version_diff.c)

$(ppdir)/version_log.c : $(sources) $(configfiles) 
#	@if hash git 2>/dev/null; then \
#		echo "Commit headlines for blessed~1..HEAD" > $(ppdir)/version_log; \
#		git log blessed~1..HEAD --pretty=format:'%h - %s (%ci)' --abbrev-commit >> $(ppdir)/version_log; \
#		echo >> $(ppdir)/version_log; \
#	else \

		sleep 10; \
		echo 1>&2; \
		echo "WARNING! No git installed on build host. Not saving any log!" 1>&2; \
		echo 1>&2; \
		echo "No git installed on build host" > $(ppdir)/version_log; \
#	fi;
	@$(call bin_to_c_file, _version_log, $(ppdir)/version_log, $(ppdir)/version_log.c)

$(ppdir)/$(notdir $(CONFIG)).c : $(configfiles)
	@$(call echo,(cat $(CONFIG)* > $(ppdir)/$(notdir $(CONFIG))))
	@$(call bin_to_c_file, _make_inc_config, $(ppdir)/$(notdir $(CONFIG)), $(ppdir)/$(notdir $(CONFIG)).c)

$(ppdir)/version_commit_id.c : $(sources) 
#	@if hash git 2>/dev/null; then \
#		git rev-parse HEAD > $(ppdir)/version_commit_id; \
#	else \

		sleep 10; \
		echo 1>&2; \
		echo "WARNING! No git installed on build host. Not saving commit id!" 1>&2; \
		echo 1>&2; \
		echo "No git installed on build host" > $(ppdir)/version_commit_id; \
#	fi;
	@$(call bin_to_c_file, _version_commit_id, $(ppdir)/version_commit_id, $(ppdir)/version_commit_id.c)

# This target will be re-made if anything changes:
$(ppdir)/compile_time.c: $(fortran_objects) $(c_objects) $(cuda_objects) $(text_objects) $(builddir)/version.o
	@date > $(ppdir)/compile_time 2>&1
	@$(call bin_to_c_file, _compile_time, $(ppdir)/compile_time, $(ppdir)/compile_time.c)

# Info in version.c
CPPFLAGS_version.c += -DMACHINE=$(HOST)
CPPFLAGS_version.c += -DUSER=$(USER)

# Python modules
pythonflags=WITH_PIC=1 WITH_OPENMP=0 WITH_HDF5=0 WITH_HDF5_SHARED=1
%.so:
	$(MAKE) --no-print-directory MAIN=$*_python $(pythonflags) $*.python

%.so.info:
	$(MAKE) --no-print-directory MAIN=$*_python $(pythonflags) info

ifeq ($(DEBUG),HIGH)
  F2PY_DEBUG = --debug
  F2PY_OPT = --noopt
endif
ifeq ($(DEBUG),MEDIUM)
  F2PY_DEBUG = --debug
  F2PY_OPT =
endif
ifeq ($(DEBUG),LOW)
  F2PY_DEBUG = --debug
  F2PY_OPT =
endif

# .python is just a make-internal tag to switch to f2py, could be .foobar just as well
%.python: $(builddir)/$(mainfile:.F90=.f90) $(filter-out %/$(mainfile:.f90=.o),$(fortran_objects)) $(c_objects) $(cuda_objects) $(text_objects) $(builddir)/compile_time.o $(builddir)/version.o
	cd build/$(id);\
	ln -sf ../../.f2py_f2cmap .; \
	../../tools/f2py.py -c $(F2PY_OPT) $(F2PY_DEBUG) \
	  --compiler=$(F2PY_CCOMPILER) \
	  --fcompiler=$(F2PY_FCOMPILER) \
	  --f90flags="$(FFLAGS)" \
	  --cflags="$(CFLAGS)" \
	  --ldflags="$(LDFLAGS)" \
	  -m $* \
	  objects/$(mainfile:.F90=.f90) \
	  -Iobjects \
	  $(patsubst $(builddir)/%, objects/%, $(filter %.o, $^))

phonies += clean
clean:
	@rm -rf build/

phonies += doc-clean
doc-clean:
	@rm -rf doc

forcheck: $(ppsources)
	@forchk -i4 -allc -decl -ancmpl -anref -shmod \
	 lib src/libs/forcheck/MPI.flb -lib src/libs/forcheck/hdf5_mpi.flb \
	 $(ppdir)/*.f90 $(ppdir)/*.f

# run a command $(P) with all preprocessed fortran
# files as argument, use as
#
#   make run-on-pp P=command
#
phonies += run-on-pp
source_of_object = $(filter $(patsubst $(builddir)/%.o, %.F90, $1), $(sources)) $(filter $(patsubst $(builddir)/%.o, %.F, $1), $(sources))
build_source = $(patsubst %.F, $(builddir)/%.f, $(patsubst %.F90, $(builddir)/%.f90, $(notdir $1)))
run-on-pp: $(foreach object,$(fortran_objects),$(call build_source, $(call source_of_object,$(object))))
	@cd $(builddir) && $(patsubst {},$(patsubst $(builddir)/%, %, $^),$(P))

phonies += doc
doc: $(call flag_dep, CPPFLAGS)
	(cat doxygen.conf; echo "INPUT = src/"; echo "PREDEFINED = $(patsubst -D%, %, $(filter -D%, $(CPPFLAGS)))") | doxygen -

# tag file for text editors
tags: $(sources)
	@files="$^";							\
	if [ "$$EDITOR" = "vi" -o "$$EDITOR" = "vim" ] ; then		\
	  ctags --globals $$files;					\
	elif [ "$$EDITOR" = "emacs" ]; then				\
	  etags --globals $$files;					\
	else								\
	  echo "Warning! Environment variable EDITOR not set!";		\
	  echo "You should always set it to your favourite editor";	\
	  echo "Generating vi-style tags file";				\
	  ctags --globals $$files;				\
	fi

# target for a single compiletest
# you can delete the output by setting
# the variable DELETE=1
phonies += $(example_configs)
$(example_configs)::
	@(echo "starting compilation test of \"$@\""; \
	  $(MAKE) --no-print-directory PREFIX=compiletest.$(@F)/ CONFIG=configurations/$@ check-cycle && \
	  $(MAKE) --no-print-directory PREFIX=compiletest.$(@F)/ CONFIG=configurations/$@ DISABLE=OPT && \
	  if [ "$(DELETE)" = 1 ] ; then \
	    rm -rf build/compiletest.$(@F)/; \
	  fi && \
	  echo -en "$@ completed without errors.\n\n") || \
	 (echo -en "\n$@ failed!\n\n"; false);

phonies += check-preprocessor
check-preprocessor: $(sources)
	@for file in $^; do \
	  tools/check_preprocessor.pl $$file || exit 1; \
	done
	@echo "No obvious preprocessor abuses"

phonies += check-cycle
check-cycle: $(depgraph)
	@$(call echo,tools/cycle_detection.pl $<)
	@echo "No cyclic dependencies"

phonies += check-no-tabs
check-no-tabs: $(filter %.F %.F90,$(sources))
	@grep -l "	" $^ && { echo TABS found in source files; exit 1; } || echo No TABS found

phonies += test
test::
	@$(MAKE) --no-print-directory check-no-tabs
	@$(MAKE) --no-print-directory check-preprocessor
	@$(MAKE) --no-print-directory check-cycle DISABLE=OPT PREFIX=compiletest.make.inc.vertex/
	@$(MAKE) --no-print-directory DISABLE=OPT PREFIX=compiletest.make.inc.vertex/
	@for compiletest in $(example_configs) remap; do \
	  echo "Testing $$compiletest" \
	  $(MAKE) --no-print-directory DELETE=0 $$compiletest || exit 1; \
	 done
	@echo -en "\nTest completed without errors.\n\n"

info_line = fmt -w 51 <<< "$(strip $($1))" | awk '{ printf("%26s%s\n", NR > 1 ? "" : "$1 = ", $$0);}' ;

phonies += defined-flags
defined-flags:
	@tr " " "\n" <<< "$(strip $(filter-out -U%,$(CPPFLAGS)))" | sort

phonies += info
info:
	@echo Architecture:
	@echo "  $(ARCH)"
	@echo
	@echo Preprocessor directives
	@$(call info_line,CPPFLAGS)
	@$(foreach file,$(sources),$(foreach mod,ADD REMOVE REPLACE,$(if $(CPPFLAGS_$(mod)_$(file)),$(call info_line,CPPFLAGS_$(mod)_$(file)),)))
	@echo
	@echo Activated compilation features:
	@echo "  $(FEATURES)"
	@echo
	@echo Compiler and link flags:
	@$(call info_line,FC)
	@$(call info_line,FFLAGS_FREE)
	@$(call info_line,FFLAGS_FIXED)
	@$(call info_line,FFLAGS)
	@$(foreach file,$(patsubst %.F90, %.f90, $(filter %.F90, $(sources))) $(patsubst %.F, %.f, $(filter %.F, $(sources))),$(foreach mod,ADD REMOVE REPLACE,$(if $(FFLAGS_$(mod)_$(file)),$(call info_line,FFLAGS_$(mod)_$(file)),)))
	@echo
	@$(call info_line,CC)
	@$(call info_line,CFLAGS)
	@$(foreach file,$(filter %.c, $(sources)),$(foreach mod,ADD REMOVE REPLACE,$(if $(CFLAGS_$(mod)_$(file)),$(call info_line,CFLAGS_$(mod)_$(file)),)))
	@echo
	@$(call info_line,CUDACC)
	@$(call info_line,CUDACFLAGS)
	@$(foreach file,$(filter %.cu, $(sources)),$(foreach mod,ADD REMOVE REPLACE,$(if $(CUDACFLAGS_$(mod)_$(file)),$(call info_line,CUDACFLAGS_$(mod)_$(file)),)))
	@echo
	@$(call info_line,CUDAFC)
	@$(call info_line,CUDAFFLAGS)
	@$(foreach file,$(filter %.cuf, $(sources)),$(foreach mod,ADD REMOVE REPLACE,$(if $(CUDAFFLAGS_$(mod)_$(file)),$(call info_line,CUDAFFLAGS_$(mod)_$(file)),)))

	@echo
	@$(call info_line,LD)
	@$(call info_line,LDFLAGS)
	@echo
	@echo Main program:
	@echo "  $(MAIN) in $(mainfile)"
	@echo
	@echo Executable:
	@echo "  $(exec)"
	@echo

phonies += why-need
why-need: $(depgraph)
	@$(if $(FILE),./tools/why_need.py $(depgraph) $(mainfile:.F90=.f90) $(FILE),echo "Usage: $(MAKE) why-need FILE=somefile.f90"; exit 1;)

phonies += executable-name
executable-name:
	@echo $(exec)

phonies += help
help:
	@grep "^##" Makefile | cut -b3-

phonies += help-features
help-features:
	@grep "^#F" Makefile | cut -b3-

phonies += list-cherries
list-cherries:
	@git fetch origin
	@red=$$(git config --get-color "" "red");			\
	green=$$(git config --get-color "" "green");			\
	blue=$$(git config --get-color "" "blue");			\
	yellow=$$(git config --get-color "" "yellow");			\
	reset=$$(git config --get-color "" "reset");			\
	list="";							\
	duplicate_before=0;						\
	branch_name="$$(git symbolic-ref HEAD 2>/dev/null)" ||		\
		branch_name="(unnamed branch)"; \
	branch_name="$${branch_name##refs/heads/}";			\
	branch_name=$${yellow}$${branch_name}$${reset};			\
	echo "You are currently on $${branch_name}";			\
	echo "The other $${yellow}user//master$${reset} branches offer the following commits not in $${branch_name}:";	\
	echo;								\
	echo "  $${red}Red commits$${reset} correspond to equivalent commits already in your branch";	\
	echo "  $${green}Green commits$${reset} are truly new revisions";	\
	echo;								\
	for b in $$(git branch -r | grep "origin/user/.*/master") ; do	\
		cherries=$$(git cherry HEAD $$b);			\
		if [ -n "$$cherries" ] ; then				\
			echo "$${green}* Branch $$b $${reset} offers the following commits:"; \
			echo "  |";					\
			while read sign id; do				\
				if [ "$${sign}" = "+" ] ; then		\
					git --no-pager show --pretty=format:'  | %Cgreen%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr)%Creset%n' --abbrev-commit -s $$id;	\
				elif [ "$${sign}" = "-" ] ; then	\
					git --no-pager show --pretty=format:'  | %Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr)%Creset%n' --abbrev-commit -s $$id;	\
				fi;					\
			done <<< "$$cherries";				\
			echo "  |";					\
			echo "  -> Consider rebasing with $${blue}git rebase $$b$${reset}";	\
			echo;						\
		fi							\
	done;

.PHONY: $(phonies)
