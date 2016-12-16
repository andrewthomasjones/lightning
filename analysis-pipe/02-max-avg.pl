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
$idir = "$ARGV[0]/00-orig";
$odir = "$ARGV[0]/";

# max of 200 frames
$nframes = 200;
@mathfiles = ();
for($t=1000; $t<(1000+$nframes); $t++){
   $ifn = sprintf("$idir/R-out-%06d.mnc", $t);

   push(@mathfiles, $ifn);
   }

&do_cmd_batch("MAX-$$-$t", 'none',
      'mincmath', '-clobber', '-nocheck_dimensions',
      '-maximum', @mathfiles,
      "$odir/02-mip.mnc");


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
