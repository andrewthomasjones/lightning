/* volcorrts.c                                                                 */
/*                                                                           */
/* A widget for correlating timeseries data                                  */
/*
/* Andrew Janke - a.janke@gmail.com                                          */
/* University of Queensland                                                  */
/*                                                                           */
/* Copyright Andrew Janke, The University of Queensland
/* Permission to use, copy, modify, and distribute this software and its     */
/* documentation for any purpose and without fee is hereby granted,          */
/* provided that the above copyright notice appear in all copies.  The       */
/* author and the University of Queensland make no representations about the */
/* suitability of this software for any purpose.  It is provided "as is"     */
/* without express or implied warranty.                                      */

/* xcorr  - Cross Correlation                                                */
/*        =  sum((a*b)^2) / (sqrt(sum(a^2)) * sqrt(sum(b^2)) */


#include <float.h>
#include <volume_io.h>
#include <ParseArgv.h>
#include <time_stamp.h>
#include <ctype.h>

#ifndef FALSE
#  define FALSE 0
#endif
#ifndef TRUE
#  define TRUE 1
#endif

#define WORLD_NDIMS 3
#define SQR2(x) ((x) * (x))

/* function prototypes */
static void print_version_info(void);

static int verbose = FALSE;
static int clobber = FALSE;
static int is_signed = FALSE;
static nc_type dtype = NC_FLOAT;

static ArgvInfo argTable[] = {
   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "General options:"},
   {"-version", ARGV_FUNC, (char *)print_version_info, (char *)NULL,
    "print version info and exit"},
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "Print out extra information."},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,
    "Clobber existing files."},

   {NULL, ARGV_HELP, NULL, NULL, "\nOutfile Options"},
   {"-byte", ARGV_CONSTANT, (char *)NC_BYTE, (char *)&dtype,
    "Write out byte data."},
   {"-short", ARGV_CONSTANT, (char *)NC_SHORT, (char *)&dtype,
    "Write out short integer data."},
   {"-long", ARGV_CONSTANT, (char *)NC_LONG, (char *)&dtype,
    "Write out long integer data."},
   {"-float", ARGV_CONSTANT, (char *)NC_FLOAT, (char *)&dtype,
    "Write out single-precision data. (Default)"},
   {"-double", ARGV_CONSTANT, (char *)NC_DOUBLE, (char *)&dtype,
    "Write out double-precision data."},
   {"-signed", ARGV_CONSTANT, (char *)TRUE, (char *)&is_signed,
    "Write signed integer data."},
   {"-unsigned", ARGV_CONSTANT, (char *)FALSE, (char *)&is_signed,
    "Write unsigned integer data."},

   {NULL, ARGV_HELP, NULL, NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
   };

char *def_4d_dimorder[] = { MIzspace, MIyspace, MIxspace, MItime };
char *def_3d_dimorder[] = { MIzspace, MIyspace, MIxspace };

int main(int argc, char *argv[]){
   char *in_fn, *out_fn;
   char *history;
   VIO_Status status;
   minc_input_options in_ops;

   VIO_Volume in_vol, out_vol;
   int in_sizes[WORLD_NDIMS], out_sizes[WORLD_NDIMS];
   double out_steps[WORLD_NDIMS], out_starts[WORLD_NDIMS];

   VIO_Real value1, value2;
   VIO_Real *t_array, *tt_array;
   int i, j, k, t, ii, jj, kk, tt;
   int c1, c2;
   double ssum1, ssum2, ssum_prd;
   double xcorr, denom;

   double max, min;

   /* get the history string */
   history = time_stamp(argc, argv);

   /* get args */
   if(ParseArgv(&argc, argv, argTable, 0) || (argc < 2)){
      fprintf(stderr,
              "\nUsage: %s [<options>] <infile.mnc> [<outfile.mnc>]\n",
              argv[0]);
      fprintf(stderr, "       %s [-help]\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }
   in_fn = argv[1];
   out_fn = argv[2];

   /* check for infile and outfiles */
   if(!file_exists(in_fn)){
      fprintf(stderr, "%s: Couldn't find input file %s.\n", argv[0], in_fn);
      exit(EXIT_FAILURE);
      }
   if(!clobber && file_exists(out_fn)){
      fprintf(stderr, "%s: File %s exists, use -clobber to overwrite.\n", argv[0],
         out_fn);
         exit(EXIT_FAILURE);
      }

   /* read in the input file */
   set_default_minc_input_options(&in_ops);
   status = input_volume(in_fn, 4, def_4d_dimorder,
      NC_UNSPECIFIED, FALSE, 0.0, 0.0, TRUE, &in_vol, &in_ops);

   if(status != VIO_OK){
      fprintf(stderr, "Problems reading: %s\n", in_fn);
      exit(EXIT_FAILURE);
      }

   if(verbose){
      fprintf(stdout, " | Input file:     %s\n", in_fn);
      }

   /* get sizes */
   get_volume_sizes(in_vol, in_sizes);

   // setup output volume
   out_sizes[0] = in_sizes[0] * in_sizes[1] * in_sizes[2];
   out_sizes[1] = in_sizes[0] * in_sizes[1] * in_sizes[2];
   out_sizes[2] = 1;

   out_steps[0] = 1;
   out_steps[1] = 1;
   out_steps[2] = 1;

   out_starts[0] = -(out_sizes[0]/2);
   out_starts[1] = -(out_sizes[1]/2);
   out_starts[2] = -(out_sizes[2]/2);

   if(verbose){
      fprintf(stdout, " | Output file:    %s\n", out_fn);
      fprintf(stdout, " | sizes  %d %d %d\n", out_sizes[0], out_sizes[1], out_sizes[2]);
      fprintf(stdout, " | steps  %g %g %g\n", out_steps[0], out_steps[1], out_steps[2]);
      fprintf(stdout, " | starts %g %g %g\n", out_starts[0], out_starts[1], out_starts[2]);
      }

   out_vol = create_volume(3, def_3d_dimorder, NC_FLOAT, TRUE, 0.0, 0.0);
   set_volume_sizes(out_vol, out_sizes);
   set_volume_separations(out_vol, out_steps);
   set_volume_starts(out_vol, out_starts);
   alloc_volume_data(out_vol);

   // for each voxel
   c1 = 0;
   min = DBL_MAX;
   max = -DBL_MAX;
   for(i = 0; i < in_sizes[0]; i++){
      for(j = 0; j < in_sizes[1]; j++){
         for(k = 0; k < in_sizes[2]; k++){


            // compare against all the other voxels
            c2 = 0;
            for(ii = 0; ii < in_sizes[0]; ii++){
               for(jj = 0; jj < in_sizes[1]; jj++){
                  for(kk = 0; kk < in_sizes[2]; kk++){


                     ssum1 = 0;
                     ssum2 = 0;
                     ssum_prd = 0;
                     for(t = in_sizes[3]; t--;){
                        GET_VALUE_4D(value1, in_vol, i, j, k, t);
                        GET_VALUE_4D(value2, in_vol, ii, jj, kk, t);

                        ssum1 += SQR2(value1);
                        ssum2 += SQR2(value2);
                        ssum_prd += SQR2(value1 * value2);
                        }

                     denom = sqrt(ssum1 * ssum2);
                     xcorr = (denom == 0.0) ? 0.0 : ssum_prd / denom;

                     set_volume_real_value(out_vol, c1, c2, 0, 0, 0, xcorr);

                     fprintf(stdout, "%g\n", xcorr);

                     if(xcorr > max){
                        max = xcorr;
                        }
                     if(xcorr < min){
                        min = xcorr;
                        }

                     c2++;
                     }
                  }
               }

            fprintf(stdout, "[%d/%d]\n", c1, out_sizes[0]);
            c1++;
            }
         }
      }

   if(verbose){
      fprintf(stdout, " | got range: %g %g\n", min, max);
      }
   set_volume_real_range(out_vol, min, max);

   if(output_volume(out_fn, dtype, is_signed, 0, 0,
                           out_vol, history, NULL) != VIO_OK){
      print_error("Problems outputing: %s", out_fn);
      }

   delete_volume(in_vol);
   delete_volume(out_vol);
   return (status);
   }

void print_version_info(void){
   fprintf(stdout, "%s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
   fprintf(stdout, "Comments to %s\n", PACKAGE_BUGREPORT);
   fprintf(stdout, "\n");
   exit(EXIT_SUCCESS);
   }
