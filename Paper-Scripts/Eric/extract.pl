#!/usr/bin/perl

# use strict;
use warnings;
#
# Eric R.Gamazon
# General extract script 

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
my $ind1 = $ARGV[2] || 0;
my $ind2 = $ARGV[3] || 0;
my %HASH1 = ();

if (!$file1 or !$file2)
{
        print "Usage: $0 <file with list of items> <file to retrieve from> <index number from first file> <index number from second file>\n";
        exit;
}

#print "Starting extraction . . .\n";

open(FH, "<$file1") or die "Can't open file1: $!";

while (my $line = <FH>)
{
        chomp($line);
        next unless ($line);

	my @data = split /\s+/, $line;
        $HASH1{$data[$ind1]} = $line;

}
close(FH);

#open(GH, "<$file2") or die "Can't open file2: $!";
open(GH, "zcat $file2 | ") or die "Can't open file2: $!"; # gzipped file

while (my $line = <GH>)
{
        chomp($line);
	$line =~ s/^\s+//;
        next unless ($line);
        my @data = split /\s+/, $line;
        if ($HASH1{$data[$ind2]})
        {
                print $line . "\n";
        }
	else
	{

#                print $line . "\t" . "0" .  "\n";
	}
}
close(GH);

