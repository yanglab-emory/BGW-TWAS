#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

Generate bfGWAS Makefile

=head1 SYNOPSIS

./gen_mkf.pl [options]

Options:
  -h      help messages
  -d      degub
  --mf    output make file

=head1 OPTIONS

=over 8

=item B<-help>

  -w         work directory : location for all output files
  -t         directory for the C++ executable file for Estep (running MCMC)
  --hyp      initial hyper parameter value file
  --em       number of EM iterations: default 5
  -b         number of burn ins: default 50,000
  -N         number of MCMC iterations: default 50,000
  --NL       number of MCMC iterations for the last EM iteration: default 50,000
  --smin     minimum number of variates per block in the model: default 0
  --smax     maximum number of variates per block in the model: default 5
  -l         specify how job will be submited (options: slurm, mosix, local): default slurm
  --mf       output make file

=back

=head1 DESCRIPTION

B<gen_mkf.pl> will generate a makefile for conducting Bayesian Functional GWASs by bfGWAS.

=cut

# define default option variables
my $help;
my $verbose;
my $debug;
my $man;
my $launchMethod = "local";
my $makeFile = "BGW-TWAS.mk";
my $genofile = "sumstat";

my $start_pos="31226778";
my $end_pos="31257725";
my $target_chr="12";
my $window_size="1e6";

my $EM=5;
my $burnin="10000";
my $Nmcmc="10000";
my $NmcmcLast="10000";

my $GTfield="GT";
my $maf="0.005";
my $rho="1";
my $smin="0";
my $smax="5";
my $win="100";
my $compress=0;
my $initype="0";
my $rv="1";
my $pp="1e-6";
my $abgamma="0.1";
# my $algorithm="VB";
my $max_iter="1000";
my $convergence="1e-7";

my $maxmem = "8000";
my $time = "24:00:00";
my $nice = "0";
my $part="nomosix";

my $seed="2020";
my $refld=0;
my $usextxld=0;
my $r2_level="0.001";
my $N="499";
my $pv="0.56";
my $BGW_dir="/home/jyang/GIT/BGW-TWAS";
my $wkDir="${BGW_dir}/Example/ExampleWorkDir";
my $hyppar="${BGW_dir}/Example/hypval.txt";
my $filelist = "${BGW_dir}/Example/ExampleData/geno_block_filehead.txt";
my $LDdir="${BGW_dir}/Example/ExampleData/LDdir";
my $Scoredir="${BGW_dir}/Example/ExampleWorkDir/ABCA7_scores";

## parameters that do not used by BGW-TWAS
my $genoDir="";
my $annoCode="";
my $pheno="";
my $annoDir="";
my $wnode="";
my $xnode="";
my $jobid="";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help, 'v'=>\$verbose, 'd'=>\$debug, 'm'=>\$man,
'w:s'=>\$wkDir, 'geno:s'=>\$genofile, 'window_size:i'=>\$window_size,
'hyp:s'=>\$hyppar, 'start:i'=>\$start_pos, 'end:i'=>\$end_pos,
'targ:s'=>\$target_chr,
'smin:s'=>\$smin, 'smax:s'=>\$smax, 'win:s'=>\$win,
'b:s'=>\$burnin, 'N:s'=>\$Nmcmc, 'NL:s'=>\$NmcmcLast,
'c:i'=>\$compress, 'n:i'=>\$N, 'converg:s'=>\$convergence, 'iter:i'=>\$max_iter,
'pp:s'=>\$pp, 'abgamma:s'=>\$abgamma,
'mem:s'=>\$maxmem, 'time:s'=>\$time, 'f:s'=>\$filelist, 'em:i'=>\$EM,
'l:s'=>\$launchMethod, 'mf:s'=>\$makeFile, 'nice:s'=>\$nice,
'j:s' =>\$jobid, 'xnode:s'=>\$xnode, 'wnode:s'=>\$wnode, 'part:s'=>\$part,
'refLD:s'=>\$refld, 'seed:s'=>\$seed, 'r2:s'=>\$r2_level,
'LDdir:s'=>\$LDdir, 'Scoredir:s'=>\$Scoredir, 'BGW_dir:s'=>\$BGW_dir,
'pv:s'=>\$pv, 'usextxLD:i'=>\$usextxld )
|| !defined($wkDir) || scalar(@ARGV)!=0)

{
    if ($help)
    {
        pod2usage(1);
        exit(0);
    }
    elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }
    else
    {
        pod2usage(1);
        exit(0);
    }
}

if ($help)
    {
        pod2usage(1);
        exit(0);
    }
elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }


if ($launchMethod ne "local" && $launchMethod ne "slurm" && $launchMethod ne "mosix" && $launchMethod ne "qsub")
{
    print STDERR "Launch method has to be local, slurm, mosix, or qsub\! \n";
    exit(1);
}

##############
#print options
##############
my $toolE="${BGW_dir}/bin/Estep_mcmc";
my $rs="${BGW_dir}/bin/Mstep.R";

printf("Options\n");
printf("\n");
printf("launch method : %s\n", $launchMethod);
printf("work directory : %s\n", $wkDir);
print "BGW-TWAS tool directory: ", $BGW_dir, "\n";
print "Estep tool: ", $toolE, "\n";
print "Mstep Rscript: ", $rs, "\n", "\n";
print "LDdir: ", $LDdir,"\n",
"Scoredir: ", $Scoredir,"\n",
"hyppar: ", $hyppar,"\n",
"fileheads: ", $filelist, "\n", "\n",
"Gene Chr: ", $target_chr, "\n",
"Gene Start Pos: ", $start_pos, "\n",
"Gene End Pos: ", $end_pos, "\n";
printf("\n");

my $comp = "";
if ($compress != 0){
  $comp = "-comp"
}

#arrays for storing targets, dependencies and commands
my @tgts = ();
my @deps = ();
my @cmds = ();

#temporary variables
my $tgt;
my $dep;
my @cmd;

mkpath($wkDir);

########################################
# Initial Set up before EM iterations
########################################

### prepare files before EM_MCMC
my $hypcurrent="$wkDir/hypval.current";
$tgt = "$wkDir/pre_em.OK";
$dep = "";
@cmd = "rm -f -r $wkDir/output $wkDir/Eoutput $wkDir/OUT";
push(@cmd, "mkdir -p $wkDir/output $wkDir/Eoutput $wkDir/OUT");
push(@cmd, "cp -f $hyppar $hypcurrent");
push(@cmd, "> $wkDir/Eoutput/EM\_result.txt");
push(@cmd, "> $wkDir/Rout.txt");
makeJob("local", $tgt, $dep, $wkDir, @cmd);


###### EM step 0 without dependencies ###########
my $i;
my $premcmcOK="";
my @filehead;
my $line;

open(my $FILELIST, $filelist)
    or die "Can not open $filelist \!";
 while ($line = <$FILELIST>) {
    chop $line;
    push(@filehead, $line);
}
close $FILELIST;

if(@filehead == 0) {
    print STDERR "file list is empty\! \n";
    exit(1);
}
else{ print "Total \# of fileheads: ", scalar(@filehead), "\n \n"; }


for(my $j=0; $j< @filehead; ++$j)
    {
        $i=0;
        $line = $filehead[$j];
        $premcmcOK .= "$wkDir/OUT/$line.$i.OK ";
        $tgt = "$wkDir/OUT/$line.$i.OK";
        $dep = "$wkDir/pre_em.OK";
        if($genofile eq "vcf"){
           @cmd = "$toolE -vcf $genoDir/$line.vcf.gz -p $pheno -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -GTfield $GTfield -maf $maf -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp  -initype $initype -rv $rv > $wkDir/OUT/$line.output.txt";
        }elsif ($genofile eq "genotxt"){
          @cmd = "$toolE -g $genoDir/$line.geno.gz -p $pheno -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -GTfield $GTfield -maf $maf -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -rv $rv > $wkDir/OUT/$line.output.txt";
        }elsif ($genofile eq "bed") {
          @cmd = "$toolE -bfile $genoDir/$line -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -GTfield $GTfield -maf $maf -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -rv $rv > $wkDir/OUT/$line.output.txt";
        }elsif($genofile eq "sumstat"){
            @cmd = "$BGW_dir/bin/run_Estep.sh $wkDir $line $N $pv $burnin $Nmcmc $LDdir $Scoredir $hypcurrent $start_pos $end_pos $target_chr $toolE ";
            # @cmd = "$toolE -score $genoDir/$line.SS.score.txt.gz -inputSS $refLD -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -n $N -pv $pv -maf $maf -r2 $r2_level -bvsrm -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -seed ${seed} > $wkDir/OUT/$line.output.txt";
        }
        else{
            die "genoDir need to be one of the vcf, genotxt, bed, or sumstat file types\!\n"
        }

        makeJob($launchMethod, $tgt, $dep, $wkDir, @cmd);
    }

my $paramfile="$wkDir/Eoutput/paramtemp$i.txt";
my $hypfile="$wkDir/Eoutput/hyptemp$i.txt";
my $logfile="$wkDir/Eoutput/log$i.txt";

$tgt = "$wkDir/Eoutput/cp_param$i.OK";
$dep = "$premcmcOK";
@cmd = "cat \`ls -d -1 $wkDir/output/** | grep paramtemp | sort\` > $paramfile";
push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep hyptemp | sort\` > $hypfile");
push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep log | sort\` > $logfile");
makeJob("local", $tgt, $dep, $wkDir, @cmd);

$tgt = "$wkDir/R$i.OK";
$dep = "$wkDir/Eoutput/cp_param$i.OK $wkDir/pre_em.OK";

@cmd = "Rscript --no-save --no-restore --verbose $rs $hypfile $i $pp $abgamma $wkDir/Eoutput/EM_result.txt $hypcurrent MCMC >> $wkDir/Rout.txt";

makeJob("local", $tgt, $dep, $wkDir, @cmd);


####### With dependencies of previous output
my $ipre="";

for $i (1..$EM){

    $ipre=$i-1; $premcmcOK="";

    for(my $j=0; $j< @filehead; ++$j){
        $line=$filehead[$j];
        $premcmcOK .= "$wkDir/OUT/$line.$i.OK ";
        $tgt = "$wkDir/OUT/$line.$i.OK";
        $dep = "$wkDir/R$ipre.OK";
        if($i < $EM){
          if($genofile eq "vcf"){
            @cmd = "$toolE -vcf $genoDir/$line.vcf.gz -a $annoDir/Anno\_$line.gz -p $pheno -fcode $annoCode -hfile $hypcurrent -GTfield $GTfield -maf $maf -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -rv $rv > $wkDir/OUT/$line.output.txt ";
          }elsif ($genofile eq "genotxt"){
            @cmd = "$toolE -g $genoDir/$line.geno.gz -p $pheno -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -GTfield $GTfield -maf $maf -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -rv $rv > $wkDir/OUT/$line.output.txt";
          }elsif ($genofile eq "bed") {
            @cmd = "$toolE -bfile $genoDir/$line -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -GTfield $GTfield -maf $maf -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -rv $rv > $wkDir/OUT/$line.output.txt";
          }elsif($genofile eq "sumstat"){
              @cmd = "$BGW_dir/bin/run_Estep.sh $wkDir $line $N $pv $burnin $Nmcmc $LDdir $Scoredir $hypcurrent $start_pos $end_pos $target_chr $toolE " ;
              # @cmd = "$toolE -score $genoDir/$line.SS.score.txt.gz -inputSS $refLD -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -n $N -pv $pv -maf $maf -r2 $r2_level -bvsrm -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -seed ${seed} > $wkDir/OUT/$line.output.txt";
          }else{
              die "genoDir need to be one of the vcf, genotxt, bed, or sumstat file types\!\n"
          }
        } elsif ($i == $EM){
            if($genofile eq "vcf"){
              @cmd = "$toolE -vcf $genoDir/$line.vcf.gz -a $annoDir/Anno\_$line.gz -p $pheno -fcode $annoCode -hfile $hypcurrent -GTfield $GTfield -maf $maf -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $NmcmcLast $comp -initype $initype -rv $rv > $wkDir/OUT/$line.output.txt ";
            }elsif ($genofile eq "genotxt"){
              @cmd = "$toolE -g $genoDir/$line.geno.gz -p $pheno -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -GTfield $GTfield -maf $maf -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $NmcmcLast $comp -initype $initype -rv $rv > $wkDir/OUT/$line.output.txt";
            }elsif($genofile eq "sumstat"){
                @cmd = "$BGW_dir/bin/run_Estep.sh $wkDir $line $N $pv $burnin $Nmcmc $LDdir $Scoredir $hypcurrent $start_pos $end_pos $target_chr $toolE ";
                #@cmd = "$toolE -score $genoDir/$line.SS.score.txt.gz -inputSS $refLD -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -n $N -pv $pv -maf $maf -r2 $r2_level -bvsrm -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -seed ${seed} > $wkDir/OUT/$line.output.txt";
            }else{
                die "genoDir need to be one of the vcf, genotxt, bed, or sumstat file types\!\n"
            }
      }
        makeJob($launchMethod, $tgt, $dep, $wkDir, @cmd);
    }

    $paramfile="$wkDir/Eoutput/paramtemp$i.txt";
    $hypfile="$wkDir/Eoutput/hyptemp$i.txt";
    $logfile="$wkDir/Eoutput/log$i.txt";

  $tgt = "$wkDir/Eoutput/cp_param$i.OK";
  $dep = "$premcmcOK";
  @cmd = "cat \`ls -d -1 $wkDir/output/** | grep paramtemp | sort \` > $paramfile";
  push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep hyptemp | sort\` > $hypfile");
  push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep log | sort\` > $logfile");
  makeJob("local", $tgt, $dep, $wkDir, @cmd);

  $tgt = "$wkDir/R$i.OK";
  $dep = "$wkDir/Eoutput/cp_param$i.OK";
  @cmd = "Rscript --no-save --no-restore --verbose $rs $hypfile $i $pp $abgamma $wkDir/Eoutput/EM_result.txt $hypcurrent MCMC >> $wkDir/Rout.txt";

    makeJob("local", $tgt, $dep, $wkDir, @cmd);

}


#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
#this tells GNU Make to delete files that were generated during a command that did not end successfully.
print MAK ".DELETE_ON_ERROR:\n\n";
#this ensures that all the targets are executed; exclude clean
print MAK "all: @tgts\n\n";

######## Create clean jobs command #######
mkpath("$wkDir/slurm_err");
push(@tgts, "clean_err");
push(@deps, "");
push(@cmds, "\t-rm -rf $wkDir/slurm_err/*.err");


push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $wkDir/*.OK $wkDir/Eoutput/*.OK $wkDir/OUT/*.OK $wkDir/slurm_err/*.err");


for(my $i=0; $i < @tgts; ++$i)
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#functions
##########

#run a job either locally or by slurm
sub makeJob
{
    my ($method, $tgt, $dep, $wkd, @cmd) = @_;

    if ($method eq "local")
    {
        makeLocalStep($tgt, $dep, @cmd);
    }
    elsif ($method eq "slurm")
    {
        makeSlurm($tgt, $dep, $wkd, @cmd);
    }
    elsif($method eq "mosix")
    {
      makeMosix($tgt, $dep, $wkd, @cmd);
    }
}

#run mosix jobs
sub makeMosix
{
    my ($tgt, $dep, $wkd, @cmd) = @_;
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmdtemp = "";
    for my $c (@cmd)
    {
      if($wnode){
        $cmdtemp .= "\tmosbatch -E$wkd -m$maxmem -j$wnode " . $c . "\n";
      }else{
        $cmdtemp .= "\tmosbatch -E$wkd -m$maxmem -b " . $c . "\n";
      }

    }
    $cmdtemp .= "\ttouch $tgt\n";
    push(@cmds, $cmdtemp);
}

#run slurm jobs
sub makeSlurm
{
    my ($tgt, $dep, $wkd, @cmd) = @_;
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmdtemp = "";
    for my $c (@cmd)
    {
        $cmdtemp .= "\tsrun --exclude=$xnode --partition=$part --mem-per-cpu\=$maxmem --time\=$time --nice\=$nice --error\=$wkDir/slurm_err/\%N.\%j.err -J $jobid -D $wkd $c \n";
    }
    $cmdtemp .= "\ttouch $tgt\n";
    push(@cmds, $cmdtemp);
}

#run a local job
sub makeLocalStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

# Check empty directories
sub dir_is_empty
{
    my ($path) = @_;
    opendir DIR, $path;
    while(my $entry = readdir DIR) {
        next if($entry =~ /^\.\.?$/);
        closedir DIR;
        return 0;
    }
    closedir DIR;
    return 1;
}
