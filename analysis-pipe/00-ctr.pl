#! /usr/bin/env perl



use warnings "all";
use File::Basename;

%opt = (
   'verbose' => 1,
   'clobber' => 1,
   'fake' => 0,
   'batch' => 1,
   'workdir' => ".",
   'tmpdir' => undef,
   );

die "Need an input dir" if !defined($ARGV[0]);
$idir = "$ARGV[0]/00-orig";


# make tmpdir
#$opt{'tmpdir'} = "./tmp/$bdir";
#&do_cmd('mkdir', '-p', $opt{'tmpdir'});
# # figure out image range
# @max = ();
# foreach $t (0,250,500,750,1000){
#
#    $t_txt = sprintf("%06d", $t);
#
#    &do_cmd('scape2mnc', '--clobber',
#       '--verbose',
#       '--floor', 0,
#       '--ceil', 50000,
#       '--start_time', $t,
#       '--stop_time', ($t+1),
#       $idir, "$opt{'tmpdir'}");
#
#    # &do_cmd('mincnorm', '-verbose',
#       # "$opt{'tmpdir'}/out-$t_txt.mnc", "$opt{'tmpdir'}/$t_txt-norm.mnc");
#
#    chomp($max = `mincstats -max -quiet $opt{'tmpdir'}/out-$t_txt.mnc`);
#    push(@max, $max);
#    }
#
# print "MAXS; " . join("|", @max) . "\n";

# get number of frames
chomp($nframes=`ls -1 $idir/out-0*.mnc | sort | tail -1`);
$nframes = &basename($nframes);
$nframes =~ s/out\-//;
$nframes =~ s/\.mnc//;
$nframes *= 1;
print "Got nframes: $nframes\n";


# centre all files
for($t=0; $t<$nframes; $t++){
   $ifname = sprintf("$idir/out-%06d.mnc", $t);
   &do_cmd_batch("VC-$$-$t", "none",
      'volcentre', $ifname);
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
