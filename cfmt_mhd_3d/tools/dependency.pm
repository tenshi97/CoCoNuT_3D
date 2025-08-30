package tools::dependency;

use strict;
use warnings;
use diagnostics;
use File::Basename;
use Data::Dumper;

# system modules
our %blacklist = (
	mpi => 1,
	ifposix => 1,
	hdf5 => 1,
	h5lt => 1,
	omp_lib => 1,
	iso_c_binding => 1
);

sub debug {
	if (defined $ENV{VERBOSE}) {
		my @args = shift;
		print STDERR @args;
	}
}

sub add_dep {
	my $dep = shift;
	my $target = shift;
	my $newdep = shift;

	if (!defined($$dep{$target})) {
		$$dep{$target} = {};
	}
	$$dep{$target}{$newdep} = 1;
}

sub compiled_name {
	my $name = shift;
	$name =~ /(.*)\.([^\.]+)/ or die;
	my $base = $1;
	my $suffix = $2;

	if ($suffix eq "F90") {
		return $base . ".f90";
	} elsif ($suffix eq "F") {
		return $base . ".f";
	} elsif ($suffix eq "f90") {
		return $base . ".o";
	} elsif ($suffix eq "f") {
		return $base . ".o";
	} elsif ($suffix eq "c") {
		return $base . ".o";
	} elsif ($suffix eq "cu") {
		return $base . ".o";
	} else {
		die "Unknown file suffix: .$suffix";
	}
}

sub scan_preprocessed_files {
	my %params = @_;
	my $ppdir = $params{ppdir};
	my $fmod = $params{fmod};

	# Map from file -> hash whose keys are modules used in the file
	my %uses;

	# Map from module -> file which provides it
	my %module_files;

	# Final result, hash of hashes
	# Map of files -> hash whose keys are prerequisites
	my %dep;

	# Find the modules used and provided by the already preprocessed files
	foreach my $file (glob($ppdir . "/{*.f90,*.f}")) {
		open(FILE,"<",$file) || die "\nCan't open $file: $!\n\n";
		$file =~ s/^$ppdir\///;
		while(<FILE>) {
			my $line = $_;
			$_ = lc($_);
			if (/^\s*module\s+procedure\b/) {
				next;
			}
			if (/^\s*use\b\s+\b([a-z0-9\-_]*)\b/) {
				add_dep(\%uses, $file, $1);
			}
			if (/^\s*module\b\s+\b(.*)\b/) {
				if (defined($module_files{$1})) {
					print STDERR "\nMODULE DEFINITION CONFLICT\n$file:$. $line\nModule '$1' also defined in file $module_files{$1}\n\n";
					exit 1;
				}
				$module_files{$1} = $file;
			}
		}
		close(FILE);
	}

	# Now that we know which file provides which module,
	# we generate the dependencies
	foreach my $file (sort keys(%uses)) {
		# Iterate over all modules the file $file needs
		foreach my $module (sort keys (%{$uses{$file}})) {
			if (!defined($module_files{$module})) {
				if (!defined($blacklist{$module})) {
					print STDERR "\nMISSING MODULE\n$file wants module '$module' which is not defined anywhere\n\n";
					exit 1;
				} else {
					next;
				}
			}
			if ($module_files{$module} ne $file) {
				add_dep(\%dep, compiled_name($file), $fmod->($module));
			}
		}
	}

	foreach my $module (sort keys %module_files) {
		add_dep(\%dep, $fmod->($module), compiled_name($module_files{$module}));
	}
	return (\%dep, \%module_files);
}

# scan unprocessed files to collect dependencies by #include
sub scan_unprocessed_files {
	my %params = @_;
	my @sources = @{$params{sources}};
	my @includedirs = sort @{$params{includedirs}};

	# Final result, hash of hashes
	# Map of files -> hash whose keys are prerequisites
	my %dep;

	# Map from program name -> file which contains it
	my %programfiles;

	# Map from program name -> List of makefile definitions
	my %programoptions;

	# Keep track of scanned files in order not to scan them twice
	my %visited;

	foreach my $file (@sources) {
		my @dependent_files = ($file);
		my %dependent_files = ($file => 1);
		foreach my $depfile (@dependent_files) {
			if (defined($visited{$depfile})) {
				debug "Skipping $depfile, already scanned\n\n";
				next;
			} else {
				$visited{$depfile} = 1;
			}
			open(FILE,"<",$depfile) || die "\nCan't open $depfile: $!\n\n";
			debug "Traversing $depfile as a dependency of $file\n";
			my @makefileoptions;
			my $program = undef;
			while(<FILE>) {
				my $line = $_;
				my $basename = $file;
				$basename =~ s/^.*\///;
				if (/^#\s*include +"(.*)"/) {
					my $included = $1;
					my $otherfile = dirname($depfile) . "/" . $included;
					# check the included file, too, if not yet done
					if(!defined($dependent_files{$otherfile})) {
						$dependent_files{$otherfile} = 1;
						push @dependent_files, $otherfile;
					}
					(my $dir = $depfile) =~ s/(.*\/).*/$1/;

					# add as dependency
					$included = $dir . $included;
					debug "$depfile:$. #include: Adding $included as a dependency of $basename\n";
					add_dep(\%dep, compiled_name($basename), $included);
				} elsif(/^#\s*include +<(.*)>/) {
					my $included = $1;

					# Only include if it is a custom include file
					# in our own paths and not a system file
					foreach my $idir (@includedirs) {
						my $otherfile = "$idir/$included";
						if (-e "$otherfile") {
							# check the included file, too, if not yet done
							if(!defined($dependent_files{$otherfile})) {
								$dependent_files{$otherfile} = 1;
								push @dependent_files, $otherfile;
							}

							# add as dependency
							add_dep(\%dep, compiled_name($basename), $otherfile);
							last;
						}
					}
				} elsif ($depfile =~ /\.(F|F90)$/ and /^\s*!make\s+(.*)$/) {
					push @makefileoptions, split(/\s+/, $1);
				} elsif ($depfile =~ /\.(F|F90)$/ and /^\s*program\b\s+\b([a-z0-9\-_]*)\b/i) {
					$program = lc($1);
					my $basename = $depfile;
					$basename =~ s/^.*\///;
					$programfiles{$program} = $basename;
				}
			}
			if (scalar(@makefileoptions) > 0) {
				die "File $depfile contained `!make' statements but no `program' section!" unless defined $program;
			}
			if (defined($program)) {
				$programoptions{$program} = \@makefileoptions;
			}
			close(FILE);
			debug "\n";
		}
	}

	return (\%dep, \%programfiles, \%programoptions);
}

1;
