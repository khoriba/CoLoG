#!/usr/bin/perl
#
# sekizuka@nih.go.jp
#
# for fastq format to one line format
#

while ($buffer = <>) {
    chomp ($buffer);
    if ($buffer =~ /^\@[MNH]/ && $buffer =~ / /)  {
    chomp ($buffer);
    print "\n".$buffer;
    }
    elsif ($buffer =~ /^\@FS/ && $buffer =~ / /)  {
    chomp ($buffer);
    print "\n".$buffer;
    }
    elsif ($buffer =~ /^\@/ && $buffer =~ / /)  {
    chomp ($buffer);
    print "\n".$buffer;
    }
    elsif ($buffer =~ /^\@[QWERTYUIOPASDFGHJKLZXCVBNM]/ && $buffer =~ / /)  {
    chomp ($buffer);
    print "\n".$buffer;
    }
    elsif ($buffer =~ /^\@read/ && $buffer =~ /\/[12]$/)  {
#    elsif ($buffer =~ /^\@[a-z0-9A-Z]/ && $buffer =~ /\/[12]$/)  {
    chomp ($buffer);
    print "\n".$buffer;
    }
    elsif ($buffer =~ /^\@/ && $buffer =~ / length/ )  {
    chomp ($buffer);
    print "\n".$buffer;
    }
#    elsif ($buffer =~ /^\@[MNHFV][0SWNFH]/)  {
    elsif ($buffer =~ /^\@[\w-_]+\:\d+\:[\w-_]+\:\d+\:\d+\:\d+\:\d+$/)  {
    chomp ($buffer);
    print "\n".$buffer;
    }
    else {
    chomp ($buffer);
    $buffer =~ s/\n//g;
    print "\t".$buffer;
    }
}



