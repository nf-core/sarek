#!/bin/sh

#Description: Create COSMIC file compatible as input to for example MuTect. Good for getting updated cosmic version.
#Requires: CosmicCodingMuts.vcf.gz and CosmicNonCodingVariants.vcf.gz from http://cancer.sanger.ac.uk/cosmic/download and the reference genome used .fai file

#Commandline: sh Create_Cosmic.sh <reference genome>.fai

#Author: Sebastian DiLorenzo

#Date: 2015-11-06

#Create the sortByRef.pl file
cat > sortByRef.pl << EOF
#!/usr/bin/perl -w

use strict;
use Getopt::Long;

sub usage {

    print "\nUsage:\n";
    print "sortByRef.pl [--k POS] [--tmp dir] INPUT REF_DICT\n\n";

    print " Sorts lines of the input file INFILE according\n";
    print " to the reference contig order specified by the\n";
    print " reference dictionary REF_DICT (.fai file).\n";
    print " The sort is stable. If -k option is not specified,\n";
    print " it is assumed that the contig name is the first\n";
    print " field in each line.\n\n";
    print "  INPUT       input file to sort. If '-' is specified, \n";
    print "              then reads from STDIN.\n";
    print "  REF_DICT    .fai file, or ANY file that has contigs, in the\n";
    print "              desired soting order, as its first column.\n";
    print "  --k POS :   contig name is in the field POS (1-based)\n";
    print "              of input lines.\n\n";
    print "  --tmp DIR : temp directory [default=/tmp]\n\n";

    exit(1);
}

my \$pos = 1;
my \$tmp = "/tmp";
GetOptions( "k:i" => \\\$pos,
	    "tmp=s" => \\\$tmp);

\$pos--;

usage() if ( scalar(@ARGV) == 0 );

if ( scalar(@ARGV) != 2 ) {
    print "Wrong number of arguments\n";
    usage();
}

my \$input_file = \$ARGV[0];
my \$dict_file = \$ARGV[1];


open(DICT, "< \$dict_file") or die("Can not open \$dict_file: \$!");

my %ref_order;

my \$n = 0;
while ( <DICT> ) {
    chomp;
    my (\$contig, \$rest) = split '\s';
    die("Dictionary file is probably corrupt: multiple instances of contig \$contig") if ( defined \$ref_order{\$contig} );

    \$ref_order{\$contig} = \$n;
    \$n++;
}

close DICT;
#we have loaded contig ordering now

my \$INPUT;
if (\$input_file eq "-" ) {
    \$INPUT = "STDIN";
} else {
    open(\$INPUT, "< \$input_file") or die("Can not open \$input_file: \$!");
}

my %temp_outputs;

while ( <\$INPUT> ) {

    my @fields = split '\s';
    die("Specified field position exceeds the number of fields:\n\$_")
        if ( \$pos >= scalar(@fields) );

    my \$contig = \$fields[\$pos];
    if ( \$contig =~ m/:/ ) {
        my @loc = split(/:/, \$contig);
        # print \$contig . " " . \$loc[0] . "\n";
        \$contig = \$loc[0]
    }
    chomp \$contig if ( \$pos == scalar(@fields) - 1 ); # if last field in line

    my \$order;
    if ( defined \$ref_order{\$contig} ) { \$order = \$ref_order{\$contig}; }
    else {
        \$ref_order{\$contig} = \$n;
        \$order = \$n; # input line has contig that was not in the dict;
        \$n++; # this contig will go at the end of the output,
              # after all known contigs
    }

    my \$fhandle;
    if ( defined \$temp_outputs{\$order} ) { \$fhandle = \$temp_outputs{\$order} }
    else {
        #print "opening \$order \$\$ \$_\n";
        open( \$fhandle, " > \$tmp/sortByRef.\$\$.\$order.tmp" ) or
            die ( "Can not open temporary file \$order: \$!");
        \$temp_outputs{\$order} = \$fhandle;
    }

    # we got the handle to the temp file that keeps all
    # lines with contig \$contig

    print \$fhandle \$_; # send current line to its corresponding temp file
}

close \$INPUT;

foreach my \$f ( values %temp_outputs ) { close \$f; }

# now collect back into single output stream:

for ( my \$i = 0 ; \$i < \$n ; \$i++ ) {
    # if we did not have any lines on contig \$i, then there's
    # no temp file and nothing to do
    next if ( ! defined \$temp_outputs{\$i} ) ;

    my \$f;
    open ( \$f, "< \$tmp/sortByRef.\$\$.\$i.tmp" );
    while ( <\$f> ) { print ; }
    close \$f;

    unlink "\$tmp/sortByRef.\$\$.\$i.tmp";
}
EOF

#Execute the Cosmic creation

#Decompress COSMIC VCFs
gunzip CosmicCodingMuts.vcf.gz
gunzip CosmicNonCodingVariants.vcf.gz

#Grab header
grep "^#" CosmicCodingMuts.vcf > VCF_Header

#Grab everything else
grep -v "^#" CosmicCodingMuts.vcf > Coding.clean
grep -v "^#" CosmicNonCodingVariants.vcf > Noncoding.clean

#Reformat to compatible input
cat Coding.clean Noncoding.clean | sort -gk 2,2 | awk '{print "chr"$0}' | perl sortByRef.pl --k 1 - $1 > Cosmic.txt

#Apply correct cosmic version information
cosmic=`grep -o -P 'COSMICv.{0,2}' VCF_Header`

#Reattach the header
cat VCF_Header Cosmic.txt > $cosmic.vcf

#clean-up directory
rm *.clean VCF_Header sortByRef.pl Cosmic.txt
