#!/usr/bin/perl 

use strict;
use warnings;

my $gmt_path = $ARGV[0];
open(my $in,  "<",  $gmt_path)  || die "Can't open input.txt: $!";

while( <$in> ){
    chomp($_);
    if ($_ =~ /^([^\t]*)\t([^\t]*)\t(.*)/){
        print "$2\t$1\t$3\n";
    } else {
        die "Error : one line doesn't match. Check your gmt.\n Line content : $_\n";
    }
}

close $in or die "$in: $!";