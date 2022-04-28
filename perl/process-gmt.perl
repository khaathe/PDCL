#!/usr/bin/perl 

use strict;
use warnings;
use Getopt::Long;

my $gmt_path;
my $min_size = 15;
my $max_size = 1000;
my $swap = 0;
my $out_path = "result.gmt";
my $help =0;

GetOptions(
    "h" => \$help,
    "min=i" => \$min_size,
    "max=i" => \$max_size,
    "swap!" => \$swap,
    "out=s" => \$out_path
) or die("Error in command line arguments\n");

if($help){
    print "Usage : ./process-gmt.perl -min <min genes> -max <max genes> -swap -out <out.gmt> input.gmt\n";
    exit(0)
}

$gmt_path = $ARGV[0];

print "Options:
    -GMT : $gmt_path
    -OUT : $out_path
    -Min : $min_size
    -Max : $max_size
    -Swap : $swap
";

$gmt_path or die("No input gmt specified");

open(my $in,  "<",  $gmt_path)  || die "Can't open $gmt_path in reading mode: $!";
open(my $out,  ">",  $out_path)  || die "Can't open $out_path in writing mode: $!";

while( <$in> ){
    chomp($_);
    if ($_ =~ /^([^\t]*)\t([^\t]*)\t(.*)/){
        my $id = $1;
        my $desc = $2;
        if ($swap){
            $id = $2;
            $desc = $1;
        }
        my $genes = $3;
        $genes =~ s/ /_/ig;
        my @gene_list = split(/\t/, $genes);
        if (@gene_list >= $min_size && @gene_list <= $max_size){
            print $out "$id\t$desc\t$genes\n";
        }
    } else {
        die "Error : one line doesn't match. Check your gmt.\n Line content : $_\n";
    }
}

close $in or die "$in: $!";
close $out or die "$out: $!";