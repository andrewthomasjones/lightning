#! /usr/bin/env perl


use warnings "all";
use File::Basename;
use File::Spec;

%opt = (
   'verbose' => 1,
   'clobber' => 1,
   'fake' => 0,
   );

die "Need an input dir" if !defined($ARGV[0]);
$me = &basename($0, (".pl"));
$input = "03-chunk02";
$idir = "$ARGV[0]/$input";
$base = &basename($ARGV[0]);
$mask = File::Spec->rel2abs("$idir/mask.nii");

&do_cmd("Rscript", "./bin/$me.R", $idir, $mask);

# run complete stats

sub do_cmd {
   print STDOUT "@_\n" if $opt{'verbose'};
   if(!$opt{'fake'}){
      system(@_) == 0 or die;
      }
   }
