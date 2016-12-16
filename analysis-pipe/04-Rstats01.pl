#! /usr/bin/env perl


use warnings "all";
use File::Basename;

%opt = (
   'verbose' => 0,
   'clobber' => 1,
   'fake' => 0,
   );

die "Need an input dir" if !defined($ARGV[0]);
$me = "04-Rstats01";
$input = "03-chunk02";
$idir = "$ARGV[0]/$input";
$odir = "$ARGV[0]";

# run R random voxel stats
$ofile = "$odir/$me.$input.pdf";
&do_cmd("Rscript", "./bin/$me.R", $idir, $ofile);


sub do_cmd {
   print STDOUT "@_\n" if $opt{'verbose'};
   if(!$opt{'fake'}){
      system(@_) == 0 or die;
      }
   }
