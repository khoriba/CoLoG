#!/usr/bin/perl
 
while ($buffer = <>) {
    if ($buffer =~ /^>/) {
    chomp ($buffer);
    $buffer =~ s/$/\%\%\%\%\%/g;
    }
    if ($buffer =~ /^\D/) {
    chomp ($buffer);
    $buffer =~ s/^>/
>/g;
    print $buffer;
    }
}
