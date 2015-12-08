#! /usr/bin/env perl



use warnings "all";
use File::Basename;

%opt = (
   'verbose' => 1,
   'clobber' => 1,
   'fake' => 0,
   'batch' => 1,
   'workdir' => ".",
   );

die "Need an input dir" if !defined($ARGV[0]);
$idir = "$ARGV[0]";
$odir = "$ARGV[0]/05-timefft";

&do_cmd('mkdir', '-p', "$odir");

# get concat files
@infiles = split(/\n/, `ls -1 $idir/????-concat.mnc`);

# get y size
chomp($ysize = `mincinfo -dimlength yspace $infiles[0]`);

# chop up
for($y=0; $y<$ysize; $y++){

   print "Doing yslice: $y\n";

   $fftfile = "$odir/" . sprintf("%03d", $y) . ".mnc";
   $fftFFT  = "$odir/" . sprintf("%03d", $y) . "-FFT.mnc";
   $fftPOW = "$odir/" . sprintf("%03d", $y) . "-POW.mnc";
   $fftMAG = "$odir/" . sprintf("%03d", $y) . "-MAG.mnc";

   @concat_files = ();
   foreach $c (@infiles){

      $ofile = "$odir/" . &basename($c, '.mnc') . "-$y.mnc";
      &do_cmd_batch("MR-$$-$y", "none",
         'mincreshape', '-clobber',
         '-dimrange', "yspace=$y,0",
         $c, $ofile);

      push(@concat_files, $ofile);
      }

   &do_cmd_batch("MC-$$-$y", "MR-$$-$y",
      'mincconcat', '-clobber', @concat_files, $fftfile);

   &do_cmd_batch("MF-$$-$y", "MC-$$-$y",
      'mincfft', '-clobber',
         "-float",
         "-1D", '-dimorder', "zspace,xspace,time",
         "-o_dimorder", "zspace,xspace,yspace",
         "-magnitude", $fftMAG,
         "-power", $fftPOW,
         $fftfile, $fftFFT);
   }


sub do_cmd {
   print STDOUT "@_\n" if $opt{'verbose'};
   if(!$opt{'fake'}){
      system(@_) == 0 or die;
      }
   }

# run a command via batch
# 1st param: Job Name
# 2nd param: Depends string
# remainder: command
# returns job ID
sub do_cmd_batch {
   my($name, $depends, $depends_str, $logdir, $buf, $jid);
   $name = shift(@_);
   $depends = shift(@_);

   $logdir = "$opt{'workdir'}/log";
   &do_cmd('mkdir', '-p', $logdir) if (!-e $logdir);

   print STDOUT "[$name:$depends] - @_\n" if $opt{'verbose'};
   if(!$opt{'fake'}){

      if($opt{'batch'}){
         print '   [B] ';
         &do_cmd('mkdir', '-p', "$logdir/$$");

         # generate and submit the script
         @args = ('qbatchf',
            '--queue', 'all.q',
            '--logfile', "$logdir/$$/$name.log",
            '--name', $name,
            '--depends', $depends,
            '--',
            @_);
         print join(' ', @args) . "\n" if $opt{'verbose'};
         &do_cmd(@args);

         print STDOUT " -- $name - $depends\n";
         }
      else{
         &do_cmd(@_);
         }
      }
   }
