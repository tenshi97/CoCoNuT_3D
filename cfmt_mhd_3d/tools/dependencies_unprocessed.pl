#!/usr/bin/perl -I./ -w

use strict;
#use warnings;
use diagnostics;

use tools::dependency;

my ($dep, $programfiles, $programoptions) = tools::dependency::scan_unprocessed_files(
	sources => [split(/ /, $ENV{SOURCES})],
	includedirs => [split(/ /, $ENV{INCLUDEDIRS})]
);

print "##############################################################\n";
print "# dependencylist.pl generated file, ". localtime(time).    " #\n";
print "##############################################################\n";
print "\n";

if (defined $ENV{VERBOSE}) {
	use Data::Dumper;
	print "# ", join("\n# ", split("\n", Dumper($programfiles))), "\n\n";
	print "# ", join("\n# ", split("\n", Dumper($dep))), "\n\n";
}

print "programs = ", join(" ", sort keys %$programfiles), "\n\n";
foreach my $program (sort keys %$programfiles) {
	print "# program ", $program, " in ", $$programfiles{$program}, "\n";
	print "mainfile_", $program, " = ", $$programfiles{$program}, "\n";
	print "\n";
	print "$program:\n";
	print "\t@\$(MAKE) --no-print-directory MAIN=$program " . join(" ", @{$$programoptions{$program}}), "\n";
	print "\n";
	print "$program-info:\n";
	print "\t@\$(MAKE) --no-print-directory MAIN=$program info " . join(" ", @{$$programoptions{$program}}), "\n";
	print "\n";
}
print "\n";

# Regex for splitting of the trailing part including the last dot
my $stem = qr/(.*)\.([^\.]+)/;

# State dependencies for all unprocessed files
foreach my $file (keys %$dep) {
	my @deps = sort keys %{$$dep{$file}};
	if (scalar(@deps) > 0) { 
		print "\$(ppdir)/$file: \\\n\t", join(" \\\n\t", @deps), "\n\n";
	}
}

print "# vi",": syntax=make\n";
