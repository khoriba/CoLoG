#!/usr/bin/perl
#
# sekizuka@nih.go.jp
#
# ×ÛPå¼³ñËXNvg (072908)
#
# FASTA`®ð
#
# >test1
# agcgtactgatcgatcgacgatcgatcgatcatcg
# gctagctagctagctacgccgcggcctaaatacgt
#
#
#@±ñÈÉ
#
#@>test1
#@a
#@g
#@c
#@g
#@t
#@a
#@c
#@t
#@g
#@a
#
# multiFASTA fileàÌòR ézñðcÉÏ·µAêUmultiFASTAðsingleFASTA`®ÅooÉµÄA
# üßÄApasteÅ¡ÉÀ×Äs­û@BsðñÉµ¼·´¶B
# zñÌ·³ªÙÈéêÍAêñmultiple alignmentµAGAPª éêÉ h - hªüÁÄ¢éæ¤ÉµÄ©çÄÑmultiFASTA`®ÉµÄ¨­ÆÇ¢©àB 
#
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





