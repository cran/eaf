#!/usr/bin/perl -w
#---------------------------------------------------------------------
#
# eafdiff.pl $Revision: 190 $
#
#---------------------------------------------------------------------
#
# Copyright (c) 2007, 2008, 2010, 2011  
# Manuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be>
# LaTeX: \copyright 2007, 2008, 2010, 2011 
#        Manuel L{\'o}pez-Ib{\'a}{\~n}ez
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You can obtain a copy of the GNU General Public License at
# http://www.gnu.org/licenses/gpl.html or writing to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
#
#---------------------------------------------------------------------
#
# IMPORTANT NOTE: Please be aware that the fact that this program is
# released as Free Software does not excuse you from scientific
# propriety, which obligates you to give appropriate credit!  If you
# write a scientific paper describing research that made substantive
# use of this program, it is your obligation as a scientist to (a)
# mention the fashion in which this software was used in the Methods
# section; (b) mention the algorithm in the References section. The
# appropriate citation is: 
#
#    Manuel López-Ibáñez, Luís Paquete, and Thomas Stützle.
#    Exploratory Analysis of Stochastic Local Search Algorithms in
#    Biobjective Optimization. In T. Bartz-Beielstein, M. Chiarandini, 
#    L. Paquete, and M. Preuss, editors, Experimental Methods for the
#    Analysis of Optimization Algorithms, pages 209–222. Springer, 
#    Berlin, Germany, 2010.
#    doi: 10.1007/978-3-642-02538-9_9
#
# Moreover, as a personal note, I would appreciate it if you would
# email <manuel.lopez-ibanez@ulb.ac.be> with citations of papers
# referencing this work so I can mention them to my funding agent and
# tenure committee.
#
#---------------------------------------------------------------------
#
# Literature:
#
# [1] Manuel Lopez-Ibanez, Luis Paquete, and Thomas Stutzle. Hybrid
#     population-based algorithms for the bi-objective quadratic
#     assignment problem.  Journal of Mathematical Modelling and
#     Algorithms, 5(1):111-137, April 2006.
#
# [2] Manuel Lopez-Ibanez, Luis Paquete, and Thomas Stutzle. Hybrid
#     population-based algorithms for the bi-objective quadratic
#     assignment problem.  Journal of Mathematical Modelling and
#     Algorithms, 5(1):111-137, April 2006.
# 
#---------------------------------------------------------------------
#
# TODO: 
#
# * Fail if any subcommand fails. For example, do not even run R if
#   the call to the eaf program failed.
#
#---------------------------------------------------------------------
use strict;
use warnings FATAL => 'all';
use diagnostics;
use File::Basename;

my $progname = basename($0);
my $version = ' $Revision: 190 $ ';

my $Rexe = "R"; 

my $ps2epsfilter = undef;
$ps2epsfilter = "ps2eps"; ## Comment this out if you do not want ps2eps.
my $eps2png = "convert"; # this is only required if --png is given

requiredprog($Rexe, $ps2epsfilter);

# These have to be kept in sync.
my $colors = '"#FFFFFF", "#BFBFBF","#808080","#404040","#000000"';
my $legend = '"[0.0, 0.2)","[0.2, 0.4)","[0.4, 0.6)","[0.6, 0.8)","[0.8, 1.0]"';

my $save_temps = 0;
my $compress_flag = 0;
my $flag_xmaximise = 0;
my $flag_ymaximise = 0;
my $fulleaf_flag = 0;
my $area_flag;
my $jpeg_flag = 0;
my $png_flag  = 0;
my $verbose_flag = 0;
my $filter = "";
my $legendpos;
my $xlim = "NULL";
my $ylim = "NULL";
my $label_left;
my $label_right;
my $label_obj1 = "objective 1";
my $label_obj2 = "objective 2";
my $output_eps;
my $overwrite = 1;

&usage(1) if (@ARGV == 0);

sub usage {
    my $exitcode = shift;
    print<<"EOF";
Usage: $progname FILE1 FILE2
Create a plot of the differences between the EAFs of FILE1 and FILE2.

 -h, --help          print this summary and exit
     --version       print version number and exit
 -v, --verbose       print some information about what is doing

     --left=STRING   label for left plot
     --right=STRING  label for right plot
     --obj1=STRING   label for objective 1 (x-axis)
     --obj2=STRING   label for objective 2 (y-axis)
                     (labels can be R expressions,
                      e.g., --obj1="expression(pi)")

     --legendpos={top,bottom}{left,right}  position of the legend

     --xlim=REAL,REAL  limits of x-axis
     --ylim=REAL,REAL  limits of y-axis

     --maximise      handle a maximisation problem
     --xmaximise     maximise first objective
     --ymaximise     maximise second objective

 -o  --output=FILE   output file
     --full          plot the full EAF instead of the differences
     --points        plot the full EAF as points instead of areas
 -z  --gzip          compress the output file (adds gz extension).
     --png           generate also png file from the eps output
                     (requires ImageMagick)
     --save-temps    Keep temporary files in the current directory
EOF
    exit $exitcode;
}

## Format commandline.
my $commandline = &format_commandline ();

my @files = ();
## Handle parameters
while (@ARGV) {
    my $argv = shift @ARGV;

    if ($argv =~ /^--left=/
           or $argv =~ /^--left$/) {
        $label_left = &get_arg ($argv);

    } elsif ($argv =~ /^--right=/
           or $argv =~ /^--right$/) {
        $label_right =  &get_arg ($argv);

    } elsif ($argv =~ /^--obj1=/
           or $argv =~ /^--obj1$/) {
        $label_obj1 = &get_arg ($argv);

    } elsif ($argv =~ /^--obj2=/
           or $argv =~ /^--obj2$/) {
        $label_obj2 =  &get_arg ($argv);

    } elsif ($argv =~ /^-o=/ or $argv =~ /^--output=/
             or $argv =~ /^-o$/ or $argv =~ /^--output$/) {
        $output_eps = &get_arg ($argv);
    }
    elsif ($argv =~ /^--legendpos=/
           or $argv =~ /^--legendpos$/) {
        my $arg = &get_arg ($argv);
        $legendpos  = "$arg" if ($arg);
    }
    elsif ($argv =~ /^--xlim=/
           or $argv =~ /^--xlim$/) {
        my $arg = &get_arg ($argv);
        $xlim  = "c($arg)" if ($arg);
    }
    elsif ($argv =~ /^--ylim=/
           or $argv =~ /^--ylim$/) {
        my $arg = &get_arg ($argv);
        $ylim  = "c($arg)" if ($arg);
    }
    elsif ($argv =~ /--maxim/) {
        $flag_xmaximise = 1;
        $flag_ymaximise = 1;
    }
    elsif ($argv =~ /--xmax/) {
        $flag_xmaximise = 1;
    }
    elsif ($argv =~ /--ymax/) {
        $flag_ymaximise = 1;
    }
    elsif ($argv =~ /--full/) {
        $fulleaf_flag = 1;
        $area_flag = 1 unless (defined($area_flag));
    }
    ## Silently accept this option for backwards compatibility.
    elsif ($argv =~ /^--area$/) {
        $area_flag = 1;
    }
    elsif ($argv =~ /^--points$/) {
        $area_flag = 0;
    }
    elsif ($argv =~ /^--png$/) {
        requiredprog($eps2png);
        $png_flag = 1;
    }
    elsif ($argv =~ /^-z$/ or $argv =~ /^--gzip$/) {
        requiredprog ("gzip");
        $compress_flag = 1;
    }
    elsif ($argv =~ /^-save-temp$/ or $argv =~ /^--save-temp/) {
        $save_temps = 1;
    }

    ## The remainder options are standard.
    elsif ($argv =~ /^-v$/ or $argv =~ /^--verbose$/) {
        $verbose_flag = "-v";
    }
    elsif ($argv =~ /^--help/ or $argv =~ /^-h/) {
        &usage(0);
    }
    elsif ($argv =~ /^--version/) {
        print "$progname: version $version\n";
        print <<'EOF';
Copyright (C) 2008-2009 Manuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be>
This is free software.  You may redistribute copies of it under the terms of
the GNU General Public License <http://www.gnu.org/licenses/gpl.html>.
There is NO WARRANTY, to the extent permitted by law.
EOF
        exit (0);
    }
    elsif ($argv =~ /^--/ or $argv =~ /^-/) {
        print "$progname: unknown option $argv\n";
        &usage(1);
    } else {
        push (@files, $argv);
    }
}

unless (@files == 2) {
    print "$progname: you must specify two input files.\n";
    print "$progname: try '$progname --help' for more information.\n";
    exit 1;
}

my $file1 = $files[0];
die "$progname: error: cannot read $file1\n" unless (-r $file1);
my $file2 = $files[1];
die "$progname: error: cannot read $file2\n" unless (-r $file2);
my $outfile  = basename($file1)."-".basename($file2);

$output_eps = (($fulleaf_flag) 
               ? "${outfile}_full.eps" : "${outfile}.eps") unless (defined($output_eps) and $output_eps);

die "$progname: error: $output_eps already exists.\n" if (-e $output_eps
                                                          and not $overwrite);

$filter = "|$ps2epsfilter -s b0 -q -l -B -O -P > " if (defined($ps2epsfilter) 
                                                       and -x $ps2epsfilter);
unless (defined($legendpos)) {
    $legendpos = ($fulleaf_flag) ? "bottomleft" : "topright";
}

$label_left = basename($file1) unless (defined($label_left) and $label_left);
$label_right= basename($file2) unless (defined($label_right) and $label_right);
$area_flag = 0 unless (defined($area_flag));


$label_left = parse_expression ($label_left);
$label_right = parse_expression ($label_right);
$label_obj1 = parse_expression ($label_obj1);
$label_obj2 = parse_expression ($label_obj2);

print "$progname: generating plot $output_eps ...\n";

my $Rfile = "$$.R";
open(R, ">$Rfile") or die "$progname: error: can't open $Rfile: $!\n";
print R <<"EOFR";
#!/usr/bin/r --vanilla
# 
# R script to plot the differences
#
# This script is executable if you have littler installed. [*]
# [*] http://code.google.com/p/littler/
#
# Input: 
#        label1    = label for first plot
#        label2    = label for second plot
#        output_eps = filename for output plot
#        legendpos = location of the legend or "" for no legend.
#
# Created by $commandline
#
# $version
#

library(eaf)

col <- c($colors)
intervals <- c($legend)
filter <- "$filter"
file.left <- "$file1"
file.right <- "$file2"
title.left  <- $label_left
title.right <- $label_right
eps.file  <- "$output_eps"
legend.pos <- "$legendpos"
maximise <- c(${flag_xmaximise}, ${flag_ymaximise})
xlab <- $label_obj1
ylab <- $label_obj2
Xlim <- $xlim
Ylim <- $ylim
full.eaf <- as.logical(${fulleaf_flag})
eaf.type <- ifelse(${area_flag}, "area", "point")

EOFR

print R <<'EOFR';

data.left <- read.data.sets (file.left)
data.right <- read.data.sets (file.right)
xlim <- range(data.left[,1], data.right[,1])
ylim <- range(data.left[,2], data.right[,2])

cat(sprintf("xlim = c(%s, %s)\n", xlim[1], xlim[2]))
cat(sprintf("ylim = c(%s, %s)\n", ylim[1], ylim[2]))
if (!is.null(Xlim)) xlim <- Xlim
if (!is.null(Ylim)) ylim <- Ylim

eps.title <- eps.file
eps.file <- paste(filter, eps.file, sep="")
postscript(eps.file, title=eps.title,
           paper = "special", horizontal=F, onefile=F,
           width=9, height=5,
           family = "Helvetica")

eafdiffplot (data.left, data.right, col = col, intervals = intervals,
              full.eaf = full.eaf, type = eaf.type, legend.pos = legend.pos,
              title.left = title.left, title.right = title.right,
              cex = 1.0, cex.lab = 1.1, cex.axis = 1.0,
              xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, maximise = maximise)
dev.off()
cat (paste("eafdiffplot:", eps.file, "\n"))

EOFR

close R;

my @args = ("$Rexe --quiet --vanilla --slave < $$.R");
system(@args) == 0
or die "$progname: error: R failed to create the plots (@args)\n";

# FIXME: EPS files tend to be huge, so sometimes pays off to convert
# them to PNG or JPEG and then back to EPS. How to fix this:
#
# Option 1: Find the proper settings of 'convert' or an alternative
# tool to do the best possible conversion. Multi-platform tools are
# better. Choose either PNG or JPEG. PNG seems the best option.
#
# Option 2: Create the output of R directly in PNG format. Then,
# conversion to EPS seems easier and pdfLaTeX handles PNG. Still, the
# proper settings for PNG would need to be investigated.
#
# Option 3: create smaller EPS files, so no conversion is
# needed. Currently, points overlapped by other points are still
# present in the EPS file and take space (and time to render). Avoid
# anything that is not shown in final plot. If the points are at
# exactly the same position, we should be able to detect that and only
# plot the point that is visible. However, a point can also be
# completely overlapped by other points in different positions (since
# points in a plot have a nonzero area). No idea how to avoid plotting
# those points. Perhaps there is some way to detec this in R itself.
#
#
if ($png_flag) {
    my $output_png = $output_eps;
    $output_png =~ s/\.eps$/.png/;
    #$output_png =~ s[/eps/][/png/];
    `$eps2png -render +antialias -density 150 -background white -flatten $output_eps $output_png`;
    #my $output_png_eps = $output_png;
    #my $output_png_eps =~ s/\.png$/_png.eps/;
    #`convert $output_png $output_png_eps`;
    #`rm -f  $output_png`;
    #`gzip --force -v $output_png_eps`;
}
#                 if($jpeg_flag) {
#                     $output_jpg = $output_eps;
#                     $output_jpg =~ s/\.eps$/.jpg/;
#                     $output_jpg =~ s[/eps/][/jpg/];
#                     $output_jpg_eps = $output_jpg;
#                     $output_jpg_eps =~ s/\.jpg$/_jpg.eps/;

#                     `convert -render +antialias -density 300 $output_eps $output_jpg`;
#                     `convert $output_jpg -compress jpeg eps2:${output_jpg_eps}`;
#                     `rm -f  $output_jpg`;
#                     `gzip --force -v $output_jpg_eps`;
#                 }

## FIXME: Do this on the fly within the R script by using pipes.
if ($compress_flag) {
    print "$progname: compressing: ";
    `gzip --force -v $output_eps`;
}

unless ($save_temps) {
    unlink("$$.R", "$$.Rout");
} else {
    print "$progname: generated R script: $$.R\n";
}



###################################
# helper sub-routines
###################################

sub parse_expression {
    my $label = shift;
    $label = "\"$label\"" if ($label !~ /^expression\(/);
    return $label;
}

sub get_arg {
    my ($option, $arg) = split (/=/, $_[0], 2);
    $arg = shift @ARGV if (not $arg);
    return $arg;
}

sub format_commandline {
    my $cmd = $0 . " ";
    for (my $i=0, my $j=25; $i < @ARGV; $i++) {
        if ($i == $j) {
            $j += 25;
            $cmd .= "\\\n# ";
        }
        if ($ARGV[$i] =~ /\s/) {
            $cmd .= " \"" . $ARGV[$i] . "\"";
        } else {
            $cmd .= " " . $ARGV[$i];
        }
    }
    return $cmd;
}

use Env '@PATH';
use Cwd;
sub requiredprog {
    my $cwd = cwd();

    foreach my $program (@_)
    {
        next if not defined $program;
        unless ($program =~ m|^.?.?/|) {
            # If no path was given
            foreach my $path (@PATH) {
                if (-e "$path/$program") {
                    $program = "$path/$program";
                    last;
                }
            }
            # Try also the current directory:
            $program = "$cwd/$program" if (-e "$cwd/$program");
        }

        die "$progname: cannot find required program `$program',"
        ." either in your PATH=\"". join(':',@PATH)
        ."\" or in the current directory`$cwd'\n"
        unless ($program =~ m|^.?.?/|);
        
        die "$progname: required program `$program' is not executable\n"
        unless (-x $program);
    }
}

sub execute {
    my $command = shift;

    if ($verbose_flag) {
        print $command ."\n";
        print `$command`;
    } else {
        `$command`;
    }
    my $exit = $? >> 8;
    die "$progname: command `$command' failed with value $exit\n" if ($exit);
}
