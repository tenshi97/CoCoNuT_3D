#!/usr/bin/perl -I./ -w

use strict;
use warnings;
use diagnostics;

use tools::dependency;

# command line argument
my $mainfile = shift;
my $ppdir = shift;
my $fmod = shift;

if ($fmod =~ /FORTRAN_MOD_FILE_NAME=(.*)/) {
	$fmod = eval $1;
} else {
	die "Invalid third argument: $fmod";
}

my ($dep, $modulefiles) = tools::dependency::scan_preprocessed_files(
	ppdir => $ppdir,
	fmod => $fmod,
);

print "##############################################################\n";
print "# dependencylist.pl generated file, ". localtime(time).    " #\n";
print "##############################################################\n";
print "\n";

if (defined $ENV{VERBOSE}) {
	use Data::Dumper;
	print "# \$dep = ", join("\n# ", split("\n", Dumper($dep))), "\n\n";
	print "# \$modulefiles = ", join("\n# ", split("\n", Dumper($modulefiles))), "\n\n";
}

# Collect all required files from the dependency graph, using the file
# with the main routine as starting point
my %reqs;
(my $object = $mainfile) =~ s/(.*)\.([^\.]+)/$1.o/;
if (defined($$dep{$object})) {
	%reqs = %{$$dep{$object}};
} else {
	%reqs = ();
}

my @reqs = keys %reqs;
foreach my $req (@reqs) {
	foreach my $subreq (keys %{$$dep{$req}}) {
		if (!defined($reqs{$subreq})) {
			push @reqs, $subreq; # recursion implemented as iteration
			$reqs{$subreq} = 1;
		}
	}
}
@reqs = sort keys %reqs;

print "fortran_objects = \$(builddir)/$object \\\n\t",
      join(" \\\n\t", map {"\$(builddir)\/$_"} grep(!/.*\.mod$/, @reqs)),
      "\n\n";

my %modules_in_file;
foreach my $module (sort keys %$modulefiles) {
	my $file = $$modulefiles{$module};
	$modules_in_file{$file} = [] unless defined $modules_in_file{$file};
	push @{$modules_in_file{$file}}, $module;
	print "\$(builddir)\/", $fmod->($module), ": \$(builddir)\/", tools::dependency::compiled_name($file), "\n";
		print "\t\@if [ \$^ -nt \$@ ] ; then echo \"Error: module file \$@ is older than object file \$^, aborting.\" 1>&2; exit 1; fi\n\n";
}

print "\n";
foreach my $file (sort keys %modules_in_file) {
	print "file_", $file, "_modules = ", join(" ", map {$fmod->($_)} @{$modules_in_file{$file}}), "\n";
}
print "\n";


# State dependencies for all required files
foreach my $file ($object, @reqs) {
	my @deps = sort keys %{$$dep{$file}};

	if (scalar(@deps) > 0) {
		print "\$(builddir)\/$file: \\\n\t", join(" \\\n\t", map{"\$(builddir)\/$_"} @deps), "\n\n";
	}
}

print "# vi",": syntax=make\n";
