#! /usr/bin/env perl



use warnings "all";
use File::Basename;

%opt = (
   'verbose' => 0,
   'clobber' => 1,
   'fake' => 0,
   'batch' => 1,
   'workdir' => ".",
   );

die "Need an input dir" if !defined($ARGV[0]);
$idir = "$ARGV[0]/02-align";
$odir = "$ARGV[0]/03-chunk";

$logdir = "$opt{'workdir'}/log";
&do_cmd('mkdir', '-p', $logdir) if (!-e $logdir);

&do_cmd('mkdir', '-p', "$odir");

# get number of frames
chomp($nframes=`ls -1 $idir/R-out-0*.mnc | sort | tail -1`);
$nframes = &basename($nframes, ".mnc");
$nframes =~ s/R\-out\-//;
$nframes *= 1;

print "GOt nframes: $nframes\n";

# mincreshape
for($t=0; $t<$nframes; $t += 1){
   $ifile = sprintf("$idir/R-out-%06d.mnc", $t);
   $ofile = sprintf("$odir/R-out-%06d.mnc", $t);

   # chunk based upon init-model.mnc
   &do_cmd_batch("chunk-$$-$t", "none",
      'mincreshape', '-clobber',
      '-dimrange', 'xspace=210,170',
      '-dimrange', 'yspace=50,70',
      '-dimrange', 'zspace=120,150',
      # '-dimrange', 'xspace=370,120',
      # '-dimrange', 'yspace=100,70',
      # '-dimrange', 'zspace=180,200',
      $ifile, $ofile);
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
