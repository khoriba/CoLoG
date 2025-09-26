#!/usr/bin/perl -w

while (<>) {
    next unless (/\t/);
    chomp;
    $i = 0;
    @list = split('\t');
    foreach $a (@list) {
        $b{$i} .= ($a . "\t");
        $i++;
    }
}
for ($j = 0; $j < $i; $j++) {
    print $b{$j} . "\n";
}
