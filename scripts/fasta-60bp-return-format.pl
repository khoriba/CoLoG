#!/usr/bin/perl


while(<>) { 
    if (/>/) {
    chomp ($_);
         print $_."\n";
    }else{
    $offset = 0;
    chomp ($_);
        while ($line = substr($_, $offset, 60)) {
        print $line, "\n";
        $offset += 60;
	}
     }
}



