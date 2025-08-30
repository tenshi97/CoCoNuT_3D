#!/usr/bin/env perl

use strict;
use warnings;
use 5.010_000;

# Replace the results of src/macros/header.h with the actual Fortran subroutines
#
# src.F90	CPP ->	pp/.f90					Perl ->	build/.f90
# -------------------------------------------------------------------------------------------
# error_if(cond, ...)	ERROR_IF_CPP("file.f90", 42)(cond, ...)		call error_if_x("file.f90:42", cond, "cond", ...)
# abort_if(cond, ...)	ABORT_IF_CPP("file.f90", 42)(cond, ...)		call abort_if_x("file.f90:42", cond, "cond", ...)
# raise_error(...)	RAISE_ERROR_CPP("file.f90, 42)(...)		call raise_abort_x("file.f90_42", ...)
# raise_abort(...)	RAISE_ABORT_CPP("file.f90, 42)(...)		call raise_abort_x("file.f90_42", ...)
# debug(x)		DEBUG_CPP("file.f90", 42)(x)			write(*,*) "file.f90:42", "x = ", x
# dump(x)		DUMP_CPP("file.f90", 42)(msd)			call dump_("x", "file.f90:42", reshape(x, (/ size(x) /)))

# regex for anything inside balanced parantheses
# i.e. something like (asd(foa()a))
my $bp = qr/(\((?>[^()]++|(?-1))*\))/;

# perl 5.8.0 version: my $bp = qr/\((?>(?>[^()]+)|(??{$bp}))*\)/;

my $cont;
my $line_length;
if ($0 =~ /post_cpp\.f90\.pl$/) {
	$cont = "&\n     &";
	$line_length = 132 - 1;
} elsif ($0 =~ /post_cpp\.f77\.pl$/) {
	$cont =  "\n     &";
	$line_length = 72;
} else {
	die "Unknown script name $0, please use this script as either post_cpp.f90.pl or post_cpp.f77.pl";
}

sub fortran {
	my $string = shift;
	die unless scalar(@_) == 0;
	my @lines;
	while (length($string) > 0) {
		my $chunk;
		if (scalar(@lines) == 0) {
			$chunk = substr $string, 0, $line_length;
		} else {
			$chunk = substr $string, 0, $line_length - length($cont) + 1;
		}
		$string = substr $string, length($chunk);
		push @lines, $chunk;
	}
	return join($cont, @lines);
}

sub cpp_file_and_line {
	my $bp = shift;
	if ($bp =~ /^"(?<file>[^"]*)", (?<linenumber>[0-9]+)$/) {
		return "\"". $+{file} . ":" . $+{linenumber} . "\"";
	} else {
		die "Malformed arguments at $ARGV:$. \"$bp\"";
	}
}

sub error_if {
	my $file_line = shift;
	my @args = @{shift(@_)};
	my $rest = shift;

	my $condition = shift(@args);
	return "call error_if_x(" . join(", ", ($file_line, $condition, "\"$condition\"", @args)) . ")$rest";
}

sub abort_if {
	my $file_line = shift;
	my @args = @{shift(@_)};
	my $rest = shift;

	my $condition = shift(@args);
	return "call abort_if_x(" . join(", ", ($file_line, $condition, "\"$condition\"", @args)) . ")$rest";
}

sub raise_error {
	my $file_line = shift;
	my @args = @{shift(@_)};
	my $rest = shift;

	return "call raise_error_x(" . join(", ", ($file_line, @args)) . ")$rest";
}

sub raise_abort {
	my $file_line = shift;
	my @args = @{shift(@_)};
	my $rest = shift;

	return "call raise_abort_x(" . join(", ", ($file_line, @args)) . ")$rest";
}

sub dump_ {
	my $file_line = shift;
	my @args = @{shift(@_)};
	my $rest = shift;

	my $array = shift @args;

	return "call dump_(" . join(", ", ("\"$array\"", $file_line, "reshape($array, (/ size($array) /))", @args)) . ")$rest";
}

sub debug_cpp {
	my $file_line = shift;
	my @args = @{shift(@_)};
	my $rest = shift;

	return "write(*,*) " . join(", ", ($file_line, "\" " . join(", ", @args) . " = \"", @args));
}

my %macros = (
	ERROR_IF_CPP	=> \&error_if,
	ABORT_IF_CPP	=> \&abort_if,
	RAISE_ERROR_CPP	=> \&raise_error,
	RAISE_ABORT_CPP	=> \&raise_abort,
	DEBUG_CPP	=> \&debug_cpp,
	DUMP_CPP	=> \&dump_,
);

my $macros = join("|", sort keys %macros);

sub process_line {
	my $line = shift;

	$line =~ s/STRING\((([^()]*($bp)?[^()]*)+)\)/"$1"/g;

	if ($line =~ /($macros)/ && $line =~ /^(?<pre>([^!]*("[^"]*")?('[^']*')?)*)(?<macro>$macros)(?<parenblocks>($bp)+)(?<rest>.*?)(?<comment>(!.*?)?)$/) {
		my $pre = $+{pre};
		my $macro = $+{macro};
		my $parenblocks = $+{parenblocks};
		my $rest = $+{rest};
		my $comment = $+{comment} . "\n";

		my @parenblocks;
		while ($parenblocks =~ /\G$bp/g) {
			push @parenblocks, substr $1, 1, -1;
		}
		die "Malformed macro" unless scalar(@parenblocks) == 2;

		my $file_line = cpp_file_and_line($parenblocks[0]);

		my @args;
		while ($parenblocks[1] =~ /\G(?<arg>($bp|[^(),]*|"[^"]*"|'[^']*')+)(?<end>,\s*|\s*$)/g) {
			push @args, $+{arg};
			last if $+{end} !~ /,/;
		}
		$line = fortran($pre . $macros{$macro}->($file_line, \@args, $rest)) . $comment;
	}

	return $line
}

while (my $_line = <>) {

	if ($_line =~ /^\s*!\s*PERL-START\s*(?<prefix>.*)/) {
		my $_code = "use strict;\n";
		$_code .= "our \$_line;\n";
		$_code .= "our \$_n;\n";
		$_code .= "our \$_line_number;\n";
		$_code .= 'local $SIG{__WARN__} = sub { die $_[0] . " at line:\n$ARGV:$_n\n$_line\n"};'."\n";
		$_code .= $+{prefix} . "\n";

		my @_lines;
		my @_ln;

		while($_line = <>) {
			if ($_line =~ /^\s*!\s*PERL-END\s*(?<suffix>.*)/) {
				$_code .= $+{suffix} . "\n return 1;";
				last;
			} elsif ($_line =~ /^\s*!\s*PERL\s*(.*)/) {
				$_code .= $1 . "\n";
			} else {
				push @_lines, process_line($_line);
				push @_ln, $.;
				$_code .= '$_line_number = '.(scalar(@_lines) - 1).'; '.
					  '$_n = $_ln['.(scalar(@_lines) - 1).']; '.
					  '($_line = $_lines['. (scalar(@_lines) - 1) . ']) =~ s/@\[\[(.*?)\]\]/eval $1/ge; print $_line;' . "\n";
			}
		}
		eval $_code or die "Unable to evalute the following perl code:\n$_code\n\nProduced:\n\n$@\n$!\n";
	} else {
		print process_line($_line);
	}
};
