#!/usr/bin/perl -w
#---------------------------------------------------------------------
#
# eafplot.pl  ($Revision: 191 $)
#
#---------------------------------------------------------------------
#
# Copyright (c) 2005, 2008, 2009, 2010
# Manuel Lopez-Ibanez  <manuel.lopez-ibanez@ulb.ac.be>
# TeX: \copyright 2005, 2008, 2009, 2010  Manuel L{\'o}pez-Ib{\'a}{\~n}ez
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
#   the call to the eaf program failed. (See execute in eafdiff.pl).
#
#---------------------------------------------------------------------
use strict;
use warnings FATAL => 'all';
use diagnostics;
use Carp;
use File::Basename;

my $progname = basename($0);
my $version = ' $Revision: 191 $ ';

my $ps2epsfilter = undef;
$ps2epsfilter = "ps2eps";
&requiredprog($ps2epsfilter);

&usage(1) if (@ARGV == 0);

sub usage {
    my $exitcode = shift;
    print <<"EOF";
$progname version $version

Usage: $progname [OPTIONS] FILES

  Plot the best/median/worst attainment surfaces for a set of input
  points.

 -h, --help          print this summary and exit;
     --version       print version number and exit;

 -v, --verbose      increase verbosity and keep intermediate files;

 -I, --iqr, --IQR     plot IQR (25%-75% percentile) instead of best and
                      worst;
 -b, --best           plot best attainment surface;
 -m, --median         plot median attainment surface;
 -w, --worst          plot worst attainment surface;
 -p, --percentiles=INT[,INT..] 
                      plot the given percentile(s) of the EAF;

   , --extra=FILE   add extra points to the plot;

   , --extra-label=STRING  label in the legend for the extra points;

     --obj1=STRING   label for objective 1 (x-axis)
     --obj2=STRING   label for objective 2 (y-axis)
                     (labels can be R expressions, 
                      e.g., --obj1="expression(pi)")

   , --legendpos={top,bottom}{left,right}  position of the legend;

     --xlim=REAL,REAL  limits of x-axis;
     --ylim=REAL,REAL  limits of y-axis;

     --area         plot the area dominated by the attainment surfaces instead
                    of lines.

EOF
    exit $exitcode;
}

## Format commandline.
my $commandline = &format_commandline ();

my @files = ();
my $percentiles;
my $verbose_flag = 1;
my $flag_ymaximise = 0;
my $flag_area = 0;
my $flag_axislog = "";
my $extra_points = "NULL";
my $extra_labels = "NULL";
my $label_obj1 = "objective 1";
my $label_obj2 = "objective 2";
my $IQR_flag = 0;
my $legendpos = "topright";
my $xlim = "NULL";
my $ylim = "NULL";
my $do_eaf = 1;

## Handle parameters
while (@ARGV) {
    my $argv = shift @ARGV;

    if ($argv =~ /^-v$/ or $argv =~ /^--verbose$/) {
        $verbose_flag = "-v";
    }
    elsif ($argv =~ /^--help/ or $argv =~ /^-h/) {
        &usage(0);
    }
    elsif ($argv =~ /^--version/) {
        print "$progname: version $version\n";
        print <<'EOF';
Copyright (C) 2010  Manuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be>
This is free software.  You may redistribute copies of it under the terms of
the GNU General Public License <http://www.gnu.org/licenses/gpl.html>.
There is NO WARRANTY, to the extent permitted by law.
EOF
        exit (0);
    }
    elsif ($argv =~ /--ymaximise/ or $argv =~ /-ymaximise/
        # For backwards compatibility
        or $argv =~ /--stoptime/ or $argv =~ /-stoptime/) {
        print "$progname: maximising y-axis\n";
        $flag_ymaximise = 1;
    }
    elsif ($argv =~ /--area/) {
        $flag_area = 1;
    }
    elsif ($argv =~ /^-b$/ or $argv =~ /^--best$/i) {
        $percentiles .= ", " if (defined($percentiles));
        $percentiles .= "0";
    }
    elsif ($argv =~ /^-m$/ or $argv =~ /^--median$/i) {
        $percentiles .= ", " if (defined($percentiles));
        $percentiles .= "50";
    }
    elsif ($argv =~ /^-w$/ or $argv =~ /^--worst$/i) {
        $percentiles .= ", " if (defined($percentiles));
        $percentiles .= "100";
    }
    elsif ($argv =~ /-I/ or $argv =~ /--iqr/i) {
        $IQR_flag = 1;
        $percentiles .= ", " if (defined($percentiles));
        $percentiles .= "25, 50, 75";
    }
    elsif ($argv =~ /^-p$/ or $argv =~ /^--percentiles=/
           or $argv =~ /^--percentiles$/) {
        $percentiles .= ", " if (defined($percentiles));
        $percentiles .= &get_arg($argv);
    }
    elsif ($argv =~ /--no-eaf/i) {
        $do_eaf = 0;
    }
    elsif ($argv =~ /^--extra=/
           or $argv =~ /^--extra$/) {
        my $arg = &get_arg ($argv);
        if ($arg and -r $arg) {
            $extra_points  = ($extra_points eq "NULL") 
            ?  "\"$arg\"" : $extra_points . ", \"$arg\"";
        } else {
            die "cannot read file $arg.\n";
        }
    }
    elsif ($argv =~ /^--extra-label=/
           or $argv =~ /^--extra-label$/) {
        my $arg = &get_arg ($argv);
        if ($arg) {
            $extra_labels  = ($extra_labels eq "NULL") 
            ?  "\"$arg\"" : $extra_labels . ", \"$arg\"";
        }
    } elsif ($argv =~ /^--obj1=/
           or $argv =~ /^--obj1$/) {
        $label_obj1 = &get_arg ($argv);

    } elsif ($argv =~ /^--obj2=/
           or $argv =~ /^--obj2$/) {
        $label_obj2 =  &get_arg ($argv);

    } elsif ($argv =~ /^--legendpos=/
           or $argv =~ /^--legendpos$/) {
        my $arg = &get_arg ($argv);
        $legendpos  = "$arg" if ($arg);
    }
    elsif ($argv =~ /^--ylog/) {
        $flag_axislog = "y";
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
    elsif ($argv =~ /^--/ or $argv =~ /^-/) {
        print "$progname: unknown option $argv\n";
        &usage(1);
    } else {
        push (@files, $argv);
    }
}

@files = &unique(@files);
die "$progname: error: no input files given\n" unless (@files);

if (defined($percentiles) and not $do_eaf) {
    die "$progname: cannot specify --no-eaf and any EAF option" .
    " (--best, --median, --worst, --iqr, --percentiles)";
}


my $filter = "";
$filter = "|$ps2epsfilter -s b0 -q -l -B -O -P > "
if (defined($ps2epsfilter) and -x $ps2epsfilter);

if ($extra_labels eq "NULL") {
    $extra_labels = $extra_points;
}

$percentiles = "NULL" unless ($do_eaf);
# Default is best, median, worst.
$percentiles = "0, 50, 100" unless (defined($percentiles));

$label_obj1 = parse_expression ($label_obj1);
$label_obj2 = parse_expression ($label_obj2);

my $Rfile = "$$.R";
open(R, ">$Rfile") or die "$progname: error: can't open $Rfile: $!\n";
print R <<"EOFR";
#!/usr/bin/r --vanilla
# 
# R script to plot attainment surfaces
#
# This script is executable if you have littler installed. [*]
# [*] http://code.google.com/p/littler/
#
# Created by $commandline
#
# $version
#

library(eaf)

filter <- "$filter"

file.extra <- list(${extra_points})
extra.legend <- c(${extra_labels})

ymaximise <- $flag_ymaximise
legend.pos <- "$legendpos"
log <- "$flag_axislog"
Xlim <- $xlim
Ylim <- $ylim
eaf.type <- ifelse(${flag_area}, "area", "point")
xlab <- $label_obj1
ylab <- $label_obj2
percentiles <- c($percentiles)
percentiles <- if (is.null(percentiles)) percentiles else sort(percentiles)
do.eaf <- as.logical($do_eaf)

EOFR

print R <<'EOFR';
if (eaf.type == "area") {
  col <- colnames
  lty <- c("solid", "dashed", "solid")
  pch <- NA
} else {
  col <- c("black", "darkgrey", "black", "grey40", "darkgrey")
  lty <- c("dashed", "solid", "solid", "solid", "dashed")
  pch <- NA
  #lty <- c("blank")
  #pch <- c(20,21,22,23,4,5)
}

epsfile <- NULL
title   <- NULL
eaffiles <- list()
data <- NULL
attsurfs <- NULL
xlim <- NULL
ylim <- NULL

EOFR

my $num = 0;
my $eaffiles = "";

foreach my $inpfile (@files) {
    unless (-r $inpfile) {
	die "$progname: warning: $inpfile: cannot read file.\n";
        #next;
    } elsif (!(-s $inpfile)
             or `grep -v -e "#\\|^\$" $inpfile | wc --bytes` =~ /\s+0$/o) {
	die "$progname: warning: $inpfile: empty file.\n";
        #next;
    }

    $num++;
    my $basefile;
    my $inpdir;
    my $file_eps;
    chomp($basefile = `basename $inpfile .dat`);
    chomp($inpdir = `dirname $inpfile`);

    my $filedata  = $inpfile;
    if ($num < 10) { print "# f$num: $inpfile\n"; }
    else           { print "#f$num: $inpfile\n"; }

    if ($do_eaf) {
        if ($IQR_flag) {
            $file_eps = "$inpdir/att_iqr_${basefile}";
        } elsif ($flag_area) {
            $file_eps = "$inpdir/att_area_${basefile}";
        } else {
            $file_eps = "$inpdir/att_${basefile}";
        }
        $eaffiles = "\"${filedata}\"";
    } else {
        my @extrafiles = `split.pl $inpfile 2>&1`;
        chomp (@extrafiles);
        $file_eps = "$inpdir/${basefile}";
        $eaffiles = "list(\"" . join ("\", \"", @extrafiles) . "\")" ;
    }
    # Include this file in the R script.
    print R <<"EOFR";

k <- $num
epsfile[k] <- paste(filter, "$file_eps", sep="")
title[k]   <- "$inpfile"
eaffiles[[k]] <- $eaffiles
EOFR

# FIXME: Use a loop in R, the loop in Perl should just fill variables.

    print R <<'EOFR';
data[[k]] <- attsurfs[[k]] <- NULL
if (do.eaf) {
  data[[k]] <- read.data.sets (eaffiles[[k]])
  xlim <- range (xlim, data[[k]][,1])
  ylim <- range (ylim, data[[k]][,2])
} else {
  attsurfs[[k]] <- lapply(eaffiles[[k]], read.data.sets)
  xlim <- range (xlim, do.call("rbind", attsurfs[[k]])[,1])
  ylim <- range (ylim, do.call("rbind", attsurfs[[k]])[,2])
}

EOFR

}

print R<<'EOFR';

if (is.null (file.extra[[1]])) {
      extra.points <-  NULL
} else {
  extra.points <- list()
  for (i in 1:length(file.extra)) {
     tmp <- read.table(file.extra[[i]], na.strings="NA")[,c(1,2)]
     extra.points[[i]] <- tmp
     xlim <- range(xlim, tmp[,1], na.rm=T)
     if (ymaximise) {
        ylim <- range(ylim, -tmp[,2], na.rm=T)
     } else {
        ylim <- range(ylim, tmp[,2], na.rm=T)
     }
  }
}

cat(sprintf("xlim = c(%s, %s)\n", xlim[1], xlim[2]))
cat(sprintf("ylim = c(%s, %s)\n", ylim[1], ylim[2]))
if (!is.null(Xlim)) xlim <- Xlim
if (!is.null(Ylim)) ylim <- Ylim

EOFR

# FIXME: Use a loop in R, the loop in Perl should just fill variables.
for (my $k = 1; $k <= $num; $k++) {
    print R<<"EOFR";

eps.file <- paste(epsfile[$k],".eps", sep='')
eps.title <- title[$k]
postscript(eps.file, title=eps.title,
           paper = "special", horizontal=F, onefile=F,
           width=4.5, height=4.5, 
           ##             width=5, height=5, 
           family = "Helvetica" ## "Times"
          )

eafplot.default (data[[$k]][,1:2], sets = data[[$k]][,3],
         attsurfs = attsurfs[[$k]], percentiles = percentiles,
         xlab = xlab, ylab = ylab, las = 0, log = log,
         type = eaf.type, lty = lty, col = col, pch=pch, cex.pch=0.75,
         extra.points=extra.points, xlim=xlim, ylim=ylim,
         legend.pos = legend.pos, extra.legend = extra.legend)

dev.null <- dev.off()
cat (paste("Plot: ",eps.file,"\\n", sep=''))

EOFR
}


close R;
&execute_verbose ("R --quiet --vanilla --slave <$Rfile");
($verbose_flag) 
? print "\n$progname: generated R script: $Rfile\n" 
: unlink($Rfile);

exit 0;

###################################
# helper sub-routines
###################################
sub parse_expression {
    my $label = shift;
    $label = "\"$label\"" if ($label !~ /^expression\(/);
    return $label;
}

sub unique {
    my %seen =() ;
    @_ = grep { ! $seen{$_}++ } @_ ;
}
sub round {
    my($number) = shift;
    return int($number + .5 * ($number <=> 0));
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

## See
## http://svn.collab.net/repos/svn/trunk/contrib/hook-scripts/svn-keyword-check.pl
## for a better way to read from a command output

sub runcmd {
    my $command = shift;
    system ($command);

    if ($? == -1) {
        die "error: failed to execute $command: $!\n";
    }
    elsif ($? & 127) {
        die "child died with signal ". ($? & 127) 
        . ($? & 128) ? "with core dump.\n" : ".\n";
    }
    else {
        return $? >> 8;
    }
}

sub execute {
    my $command = shift;

    if ($verbose_flag) {
        &execute_verbose ("$command");
    } else {
        `$command`;
    }
}

sub execute_verbose {
    unless (@_) {
        croak "$progname: execute_verbose passed no arguments.\n";
    }
    print "\n@_\n";
    my $fh = _pipe(@_);
    my @output;
    while (<$fh>) {
        print;
        chomp;
        push(@output, $_);
    }
    close($fh);
    my $result = $?;
    my $exit   = $result >> 8;
    my $signal = $result & 127;
    my $cd     = $result & 128 ? "with core dump" : "";
    if ($signal or $cd) {
        warn "$progname: pipe from `@_' failed $cd: exit=$exit signal=$signal\n";
    }
    if (wantarray) {
        return ($result, @output);
    } else {
        return $result;
    }
}

# Return the filehandle as a glob so we can loop over it elsewhere.
sub _pipe {
    local *SAFE_READ;
    my $pid = open(SAFE_READ, '-|');
    unless (defined $pid) {
        die "$progname: cannot fork: $!\n";
    }
    unless ($pid) {
        open(STDERR, ">&STDOUT") or die "$progname: cannot dup STDOUT: $!\n";
        exec(@_) or die "$progname: cannot exec `@_': $!\n";
    }
    return *SAFE_READ;
}

