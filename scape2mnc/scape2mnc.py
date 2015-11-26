#! /usr/bin/env python
#
# convert SCAPE data to a series of MINC files
#
#
# Andrew Janke - a.janke@gmail.com
#
# Copyright Andrew Janke, The University of Queensland.
# Permission to use, copy, modify, and distribute this software and its
# documentation for any purpose and without fee is hereby granted,
# provided that the above copyright notice appear in all copies.  The
# author and the University make no representations about the
# suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
#
# sudo apt-get install python-scipy

import sys
import argparse
import scipy.io as sio
import math
import os.path
import datetime
import codecs   # darn UTF-8 BOM windows encodings...
import glob
import os
import tempfile
import subprocess
from functools import partial
import time
import pprint

BLOCKSIZE = 1048576


# sort based upon file naming of ANDOR data
# eg: 7000000000spool.dat the digits are reversed, there
# is no doubt a good reason for this.
def scape_sort(item):
   # remove bits
   item.replace('spool.dat','')

   # return the reverse string
   return item[::-1]

def do_cmd(cmd):
   if args.verbose:
       print " ".join(str(x) for x in cmd)
   if not args.fake:
      subprocess.call(str(x) for x in cmd)

def mcomplete(minc_file):
    FNULL = open(os.devnull, 'w')
    status = subprocess.call(["minccomplete", minc_file], stdout=FNULL, stderr=subprocess.STDOUT)

    if(status == 0):
        return 1
    else:
        return 0

def read_blocks(files):
    for filename in files:
        if args.verbose:
           print " + reading " + filename
        with open(filename, "rb") as f:
            for block in iter(partial(f.read, BLOCKSIZE), b''):
                yield block


me = os.path.basename(sys.argv[0])

# get history string
history = datetime.datetime.now().isoformat()
history += '>>>> ' + " ".join(str(x) for x in sys.argv)

# set up and parse command line arguments
parser = argparse.ArgumentParser(description = "read a SCAPE directory and " +
                    "output to a series of MINC volumes, one per timeslice")
parser.add_argument('-v', '--verbose', help="be verbose",
                    action="store_true", default=True)
parser.add_argument('-c', '--clobber', help="clobber existing files",
                    action="store_true", default=False)
parser.add_argument('-f', '--fake', help="do a dry run (echo cmds only)",
                    action="store_true", default=False)
parser.add_argument('--floor', type=float, help="lower bound of input data")
parser.add_argument('--ceil', type=float, help="upper bound of input data")
parser.add_argument('--start_time', type=int, help="time slice to start at")
parser.add_argument('--stop_time', type=int, help="time slice to stop at")
parser.add_argument("indir", help="the input SCAPE directory")
parser.add_argument("outdir", help="the output MINC/HDF5 directory")
args = parser.parse_args()

print "Indir is " + args.indir
print "History is " + history

# make output directory if it doesn't exist
try:
    os.makedirs(args.outdir)
except:
    print args.outdir + " exists, using that"

# make tmpdir
tmpdir = tempfile.gettempdir()

# load MATLAB parameters
mat_contents = sio.loadmat(args.indir + "_info.mat")
framerate = mat_contents['info']['camera'][0][0]['framerate'][0][0][0][0]
scanrate = mat_contents['info']['daq'][0][0]['scanRate'][0][0][0][0]

if args.verbose:
    print "MAT contents"
    pprint.pprint(mat_contents)

# parse the .ini parameters
ini_data = {}
ini_meta = args.indir + "/acquisitionmetadata.ini"
section = "none"
for line in open(ini_meta, "rU"):
    if line.startswith(codecs.BOM_UTF8):
        line = line[3:]

    # strip unintersting(blank) lines
    line = line.strip().rstrip()

    # skip blank lines
    if not line:
        continue

    if line[0] == '[':
        section = line[1:-1]
    else:
        name, val = line.split(' = ')
        print "[" + section + "] NAME: " + name + " VAL: " + str(val)
        ini_data[section, name] = val

if args.verbose:
    print "INI contents"
    pprint.pprint(ini_data)

# setup parameters for conversion
xsize = int(ini_data['data', 'AOIWidth'])
ysize = int(ini_data['data', 'AOIHeight']) + 2
zsize = int(math.floor(framerate/scanrate))
tsize = 4800

# time slices to get
if args.start_time is None:
    start_time = 0
else:
    start_time = args.start_time

if args.stop_time is None:
    stop_time = tsize
else:
    stop_time = args.stop_time

if args.verbose:
   print " + xsize: " + str(xsize) + " ysize: " + str(ysize) + " zsize: " + str(zsize)
   print " + time range: [" + str(start_time) + ":" + str(stop_time) + "]"

print "Pixel Encoding " + ini_data['data', 'PixelEncoding']
if ini_data['data', 'PixelEncoding'] == "Mono16":
   convert_datatype = ('-unsigned', '-short')
else:
   raise Exception( me + ": Unknown datatype: " + ini_data['data', 'PixelEncoding'])


# figure the number of bytes per volume
bytes_per_voxel = 2
bytes_per_volume = xsize * ysize * zsize * bytes_per_voxel

# get the list of files
spoolfiles = glob.glob(args.indir + "/*.dat")

# sort them
spoolfiles.sort(key=scape_sort)
if args.verbose:
    print "Found " + str(len(spoolfiles)) + " files"
    # print spoolfiles

# get size of the first (assume for the others)
bytes_per_spool = os.path.getsize(spoolfiles[0])
if args.verbose:
    print "First spool " + spoolfiles[0] + " size is " + str(bytes_per_spool)


# create the MINC files, one volume per time step
starts = (0,0,0,0);
steps = (1,1,1,0.2);
for t in range(start_time, stop_time):
    # check how long this will take
    elapsed_start = time.time()

    outfile = args.outdir + "/out-%06d.mnc" % t

    # figure out the offset, start, stop and seek
    offset = t * bytes_per_volume
    start_spool_idx = int(math.floor(offset / bytes_per_spool))
    seek = offset % bytes_per_spool
    stop_spool_idx = int(start_spool_idx + math.floor((bytes_per_volume + seek) / bytes_per_spool)) + 1

    if args.verbose:
        print "[" + str(t) + "/" + str(tsize) + "] Offset in bytes: " + str(offset) + " + " + str(bytes_per_volume)
        print " - Using spools [" + str(start_spool_idx) + ":" + str(stop_spool_idx) + "] " + \
            str(len(spoolfiles)) + " spoolfiles @ " + str(bytes_per_spool) + "bytes"
        print " - Seek in first: " + str(seek)

    # check if we are stepping off the end of the world
    if start_spool_idx > len(spoolfiles):
        print "XXXXXXXX probably done, beyond the end of space " + str(start_spool_idx) + " > " + str(len(spoolfiles))


    # check if already done
    if(mcomplete(outfile)):
        print " - exists, skipping"

    else:
        cmd = ['rawtominc', '-clobber',
            '-zyx',
            '-xstep', str(steps[0]), '-ystep', str(steps[1]), '-zstep', str(steps[2]),
            '-xstart', str(starts[0]), '-ystart', str(starts[1]), '-zstart', str(starts[2]),
            '-dattribute', "time:start=" + str(t * steps[3]),
            '-dattribute', "time:step=" + str(steps[3]),
            '-unsigned', '-short',
            '-range', str(args.floor), str(args.ceil),
            '-orange', str(args.floor), str(args.ceil),
            '-skip', str(seek),
            outfile,
            str(1), str(zsize), str(ysize), str(xsize)]

        #print cmd
        if args.verbose:
            print "CMD: " + " ".join(str(x) for x in cmd)

        raw = subprocess.Popen(cmd, stdin=subprocess.PIPE)
        for block in read_blocks(spoolfiles[start_spool_idx:stop_spool_idx]):
            try:
                raw.stdin.write(block)
            except:
                print " - rawtominc appears full " + str(stop_spool_idx)

    elapsed_time = time.time() - elapsed_start
    print " = Done -- took " + str(elapsed_time) + "s " + \
        str((stop_time - t - 1) * elapsed_time /60) + "m to complete"
