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
$bdir = "$ARGV[0]";
$idir = "$bdir/00-orig";
$odir = "$bdir/02-align";
$model = "init-model.mnc";

$logdir = "$opt{'workdir'}/log";
&do_cmd('mkdir', '-p', $logdir) if (!-e $logdir);

&do_cmd('mkdir', '-p', "$odir");

# flip z + 2z
&do_cmd('param2xfm', '-clobber',
   '-scales', 1, 1, -1,
   '-translation', 0, 0, 2,
   "$odir/fz+2z.xfm"
   );

# shear
&do_cmd('param2xfm', '-clobber',
   '-shears', 0, 0, -1,
   "$odir/sh-1z.xfm"
   );

# flip x
&do_cmd('param2xfm', '-clobber',
   '-scales', -1, 1, 1,
   "$odir/fx.xfm"
   );

&do_cmd('xfmconcat', '-clobber',
   "$odir/sh-1z.xfm", "$odir/fx.xfm",
   "$odir/even-init.xfm");

# create MIP
if(&mcomplete("$bdir/02-align.mip.mnc")){
   printf " + $bdir/02-align.mip.mnc Exists, skipping\n";
   }
else{

   $nframes = 200;
   @mathfiles = ();
   for($t=1000; $t<(1000+$nframes); $t+=2){
      $ifn = sprintf("$idir/out-%06d.mnc", $t);
      push(@mathfiles, $ifn);
      }

   &do_cmd('mincmath', '-clobber',
      '-nocheck_dimensions',
      '-maximum', @mathfiles,
      "$odir/02-mip.mnc");
   &do_cmd('mincreshape', '-clobber',
      '-dimrange', 'time=0,0',
      "$odir/02-mip.mnc", "$odir/02-mipt.mnc");
   &do_cmd('mincresample', '-clobber',
      '-tfm_input_sampling',
      '-transformation', "$odir/sh-1z.xfm",
      "$odir/02-mipt.mnc", "$odir/02-mipt.res.mnc");
   chomp($extents = `volextents $odir/02-mipt.res.mnc`);
   &do_cmd('mincresample', '-clobber',
      '-use_input_sampling',
      split(/\ /, $extents),
      '-xdircos', 1, 0, 0,
      '-ydircos', 0, 1, 0,
      '-zdircos', 0, 0, 1,
      '-transformation', "$odir/even-init.xfm",
      "$odir/02-mipt.mnc", "$bdir/02-align.mip.mnc");
   }


# create mask
if(&mcomplete("$bdir/02-align.mask.mnc")){
   printf " + $bdir/02-align.mask.mnc Exists, skipping\n";
   }
else{
   $nframes = 200;
   @mathfiles = ();
   for($t=1000; $t<(1000+$nframes); $t+=2){
      $ifn = sprintf("$idir/out-%06d.mnc", $t);
      push(@mathfiles, $ifn);
      }

   # first create an average
   &do_cmd('mincaverage', '-clobber',
      '-nocheck_dimensions',
      @mathfiles,
      "$odir/02-avg.mnc");
   &do_cmd('mincreshape', '-clobber', '-dimrange', 'time=0,0',
         "$odir/02-avg.mnc", "$odir/02-avgt.mnc");

   # get threshold
   &do_cmd('mincnorm', '-clobber',
      "$odir/02-avgt.mnc", "$bdir/02-align.avg.mnc");

   #chomp($thresh = `mincstats -quiet -biModalT $odir/02-avgt.norm.mnc`);
   #print "Got threshold of $thresh\n";
   #if($thresh > 40){
   #   $thresh = 12;
   #}
   die "No threshold file ($bdir/00-orig.thresh.txt) found\n" .
      "   Create using -- register $bdir/02-align.avg.mnc\n"
      if !-e "$bdir/00-orig.thresh.txt";
   chomp($thresh = `cat $bdir/00-orig.thresh.txt`);
   print "Using treshold of $thresh\n";

   # reshape to MIP data
   &do_cmd('mincresample', '-clobber',
         '-like', "$bdir/02-align.mip.mnc",
         '-transformation', "$odir/even-init.xfm",
         "$bdir/02-align.avg.mnc", "$odir/02-avgt.res-sq.mnc");

   # pad, binarise, erode, group, dilate, blur, un-pad
   &do_cmd('volpad', '-clobber',
      '-distance', 10,
      "$odir/02-avgt.res-sq.mnc", "$odir/02-avgt.res-sq.pad.mnc");
   &do_cmd('mincmorph', '-clobber',
      '-successive', "B[$thresh:100:1:0]EGK[0:1]DDDDD",
      "$odir/02-avgt.res-sq.pad.mnc", "$odir/02-avgt.res-sq.morph.mnc");
   &do_cmd('mincblur', '-clobber',
      '-fwhm', 10,
      "$odir/02-avgt.res-sq.morph.mnc", "$odir/02-avgt-b");
   &do_cmd('volpad', '-clobber',
      '-distance', -10,
      "$odir/02-avgt-b_blur.mnc", "$bdir/02-align.mask.mnc");

   # check image
   &do_cmd('mph', "$bdir/02-align.mask.mnc");
   }


# align to self
if(&mcomplete("$bdir/align.xfm")){
   print " + $bdir/align.xfm Exists, skipping\n";
   }
else{
   # mask data
   &do_cmd('mincmath', '-clobber',
      '-mult', "$bdir/02-align.mip.mnc", "$bdir/02-align.mask.mnc",
      "$odir/02-reg-source.mnc");

   # initial rotation
   &do_cmd('param2xfm', '-clobber',
      '-rotations', -45, 0, 0,
      "$odir/rot--45x.xfm"
      );

   &do_cmd('mincresample', '-clobber',
      '-use_input_sampling',
      '-transformation', "$odir/rot--45x.xfm",
      "$odir/02-reg-source.mnc", "$odir/02-reg-source-rot.mnc");

   # align data
   &do_cmd('volalign', '-clobber',
      '-y',
      "$odir/02-reg-source-rot.mnc", "$odir/align.xfm", "$bdir/02-align.align.mnc");

   # add initial rotation back on
   &do_cmd('xfmconcat', '-clobber',
      "$odir/rot--45x.xfm", "$odir/align.xfm", "$bdir/align.xfm");

   # check image
   &do_cmd('mph', "$bdir/02-align.align.mnc");
   }


# register to model
$lin_xfm = "$bdir/model.lin.xfm";
if(&mcomplete($lin_xfm)){
   print " + $lin_xfm Exists, skipping\n";
   }
else{
   # register data
   &do_cmd('minctracc', '-clobber', '-debug',
      '-lsq12',
      '-transformation', "$bdir/align.xfm",
      "$odir/02-reg-source.mnc", $model, "$bdir/model.lin.xfm");

   # resample lin fit
   &do_cmd('mincresample', '-clobber',
      '-transformation', "$bdir/model.lin.xfm",
      '-like', $model,
      "$odir/02-reg-source.mnc", "$bdir/model.lin.mnc");

   &do_cmd('mph', "$bdir/model.lin.mnc");

   # resample mask
   &do_cmd('mincresample', '-clobber',
      '-transformation', "$bdir/model.lin.xfm",
      '-like', $model,
      "$bdir/02-align.mask.mnc", "$bdir/model-mask.lin.mnc");

   # check image
   &do_cmd('mph', "$bdir/model-mask.lin.mnc");
   }



# create even and odd xfms
&do_cmd('xfmconcat', '-clobber',
   "$odir/sh-1z.xfm", "$odir/fx.xfm", "$bdir/model.lin.xfm",
   "$odir/even.xfm");
&do_cmd('xfmconcat', '-clobber',
   "$odir/fz+2z.xfm", "$odir/sh-1z.xfm", "$odir/fx.xfm", "$bdir/model.lin.xfm",
   "$odir/odd.xfm");


# get number of frames
chomp($nframes=`ls -1 $idir/out-0*.mnc | sort | tail -1`);
$nframes = &basename($nframes);
$nframes =~ s/out\-//;
$nframes =~ s/\.mnc//;
$nframes *= 1;

print "GOt nframes: $nframes\n";


# re-align
for($t=0; $t<$nframes; $t+=2){
   $ifname0 = sprintf("$idir/out-%06d.mnc", $t);
   $ifname1 = sprintf("$idir/out-%06d.mnc", ($t+1));
   $ofname0 = sprintf("$odir/R-out-%06d.mnc", $t);
   $ofname1 = sprintf("$odir/R-out-%06d.mnc", ($t+1));

   if(&mcomplete($ofname0)){
      print " + $ofname0 Exists, skipping\n";
      }
   else{
      &do_cmd_batch("RES-$$-$t", "none",
         'mincresample', '-clobber',
         '-like', $model,
         '-transformation', "$odir/even.xfm",
         $ifname0, $ofname0);
      }

   if(&mcomplete($ofname1)){
      print " + $ofname1 Exists, skipping\n";
      }
   else{
      &do_cmd_batch("RES-$$-" . ($t+1), "none",
         'mincresample', '-clobber',
         '-like', $model,
         '-transformation', "$odir/odd.xfm",
         $ifname1, $ofname1);
      }
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
   my($name, $depends, $depends_str, $buf, $jid);
   $name = shift(@_);
   $depends = shift(@_);

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

# little function to test if a minc or xfm file is complete (and exists)
sub mcomplete {
   my $infile = shift(@_);

   chomp(my $buf = `minccomplete -error_string 1 $infile`);

   return ($buf == 0) ? 1 : 0;
   }
