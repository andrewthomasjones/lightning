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
$input = "04-Splines";
$input2 = "05-Clusters";
$output = "06-PostProc";
$idir = "$ARGV[0]/$input";
$idir2 = "$ARGV[0]/$input2";
$fold = (split "/", $ARGV[0])[4]; 
$odir = "/data/nif02/uqajon14/$fold/$output";
$base = &basename($ARGV[0]);
$mask = File::Spec->rel2abs("$idir/mask.nii");


# run R random voxel stats
chomp($from=`grep $base /data/nif02/uqajon14/00_FROM_TO.txt | cut -f2 -d' ' `);
chomp($to=`grep $base /data/nif02/uqajon14/00_FROM_TO.txt | cut -f3 -d' ' `);
&do_cmd("Rscript", "./temp/$me.R", $idir, $odir, $from, $to, $mask, $idir2);

# run complete stats

sub do_cmd {
   print STDOUT "@_\n" if $opt{'verbose'};
   if(!$opt{'fake'}){
      system(@_) == 0 or die;
      }
   }
