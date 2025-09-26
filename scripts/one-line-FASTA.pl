#!/usr/bin/perl
#
# new version! FASTAのフォーマットを一発で一行にするスクリプト。
#  
#　こんなやつが、
# >E8Y3Q8Q03FMQXE
# GTGATATATGTGGTGTGTGATGTGTAAGGTGTGTATGTG
# TGTAATGTGTGCATAGTGTTTGTGGTGCATGTGTGTGTG
# TGGTGTGT
#
# こんな風に！
# >E8Y3Q8Q03FMQXE:GTGATATATGTGGTGTGTGATGTGTAAGGTGTGTATGTGTGTAATGTGTGCATAGTGTTTGTGGTGCATGTGTGTGTGTGGTGTG
#
#
# nucleotide, amino acid両方のシークエンス配列に使用出来ます。
#
#
# sekizuka@nih.go.jp
#
#



 
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
