#!/usr/bin/perl
# % ./fasta_to_vertical-type.pl 'input-data' > 'output-data-1'
# % ./multi_to_single-FASTA.pl 'input-data (output-data-1)' 
# % paste 'input-data (all output-data from ./multi_to_single-FASTA.pl)'
#


while ($buffer = <>) {
    if ($buffer =~ /^>/)  {
    chomp ($buffer);
    print $buffer."
";
    } else {
    chomp ($buffer);
    @seq = split (//, $buffer);
    $base = join('
', @seq);
    print $base."
";
    }
}





