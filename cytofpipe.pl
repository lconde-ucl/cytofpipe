#!/usr/bin/env perl

use strict;
use warnings;

no warnings qw/uninitialized/;

use Getopt::Long qw(GetOptionsFromArray :config pass_through no_auto_abbrev bundling);
use File::Basename;


#######################################################################
##			ENV variables				     ##
##   Modify CYTOFPIPE_HOME to point to cytofpipe master directory    ##
#######################################################################

$ENV{'CYTOFPIPE_HOME'} = '/home/regmond/Scratch/cytof/standalone_cytofpipe/cytofpipe_v1.0';

my $os=$^O;
if($os eq 'darwin'){
	$ENV{'R_LIBS'} = $ENV{'CYTOFPIPE_HOME'}."/Rlibs_mac";
}elsif ($os eq 'linux'){
	$ENV{'R_LIBS'} = $ENV{'CYTOFPIPE_HOME'}."/Rlibs";	
}else{
	$ENV{'R_LIBS'} = $ENV{'CYTOFPIPE_HOME'}."/Rlibs";	
}
$ENV{'R_MAX_NUM_DLLS'} = 153;
$ENV{'RAND_ID'} = `od -N 4 -t uL -An /dev/urandom | tr -d " " | tr -d "\n"`;
#######################################################################


my @arg0=($ARGV[0]);

my $command='';
GetOptionsFromArray (
    \@arg0,
    "clustering" => \&parse_clustering,
    "scaffold" => \&parse_scaffold,
    "citrus" => \&parse_citrus,
    "<>"   => \&print_usage
) or die "\n";


sub parse_clustering {

	my @args=@ARGV;
	my $command=shift(@args);

	my $inputdir=''; my $outputdir=''; my $markersfile='';
	my $configfile='';
	my $flow='0';
	my $cytof='0';
	my $displayall='';
	my $all='';
	my $downsample='';

	GetOptionsFromArray (
	    \@args,
	    "i=s" => \$inputdir,
	    "o=s"   => \$outputdir,
	    "m=s"   => \$markersfile,
	    "config=s"   => \$configfile,
	    "flow"   => \$flow,
	    "cytof"   => \$cytof,
	    "displayAll"   => \$displayall,
	    "all"   => \$all,
	    "downsample=i"   => \$downsample,
	    "<>"   => \&print_clustering
	) or die "\n";

        if ($inputdir eq '' || $outputdir eq '' || $markersfile eq ''){
                usage_clustering("Please check that you are providing a inputdir (-i), outputdir (-o) and markersfile (-m)");
		return;
        }
	if (!-e "$inputdir") {
                usage_clustering("Can't find directory with fcs files <$inputdir>");
		return;
	}		
	if (-e "$outputdir") {
		usage_clustering("The outputdir <$outputdir> already exists, please choose a different outputdir");
		return;
	}
	if (!-e "$markersfile") {
		usage_clustering("Can't find markers file <$markersfile>");
		return;
	}
	if($flow eq '1' && $cytof eq '1') {
		usage_clustering("These two parameters [--flow, --cytof] can not be used jointly, please choose one of them (or none for default options)");
		return;
	}
	if($downsample ne ''){
		if($all eq '1' && $downsample ne ''){
			usage_clustering("These two parameters [--all, --downsample NUM] can not be used jointly, please choose one of them (or none for default options)");
        	      	return;
		}
		if(!isnum($downsample) || ($downsample < 500 || $downsample > 100000)){
			usage_clustering("<$downsample> is not a valid downsample option. Please insert an integer between 500 and 100,000");
			return;
		}	
	}
	if ($configfile ne ''){
		if (!-e "$configfile") {
			usage_clustering("Can't find config file <$configfile>");
			return;
		}else{
			check_config_clustering($configfile)
		}
	}

	my $run="$ENV{CYTOFPIPE_HOME}/code/run_cytofpipe.sh "."@ARGV";
	system("$run");
}



sub parse_scaffold {

	my @args=@ARGV;
        my $command=shift(@args);


	my $inputdir=''; my $outputdir=''; my $markersfile='';
	my $reference="";
	my $flow=0;
	my $cytof=0;
	my $all='';
	my $downsample='';
	GetOptionsFromArray (
	    \@args,
	    "i=s" => \$inputdir,
	    "o=s"   => \$outputdir,
	    "m=s"   => \$markersfile,
	    "ref=s"   => \$reference,
	    "flow"   => \$flow,
	    "cytof"   => \$cytof,
	    "all"   => \$all,
	    "downsample=i"   => \$downsample,
	    "<>"   => \&print_scaffold
	) or die "\n";

	if ($inputdir eq '' || $outputdir eq '' || $markersfile eq '' || $reference eq '') {
		usage_scaffold("Please check that you are providing a inputdir (-i), outputdir (-o), markersfile (-m) and reference FCS (--ref)");
		return;
	}		

	if (!-e "$inputdir") {
		usage_scaffold("Can't find directory with fcs files <$inputdir>");
		return;
	}		
	if (!-e "${inputdir}/gated") {
		usage_scaffold("Can't find the \'gated/\' subfolder inside inputdir <$inputdir>");
	       	return;
	}
	if (!-e "${inputdir}/${reference}") {
		usage_scaffold("Can't find the reference FCS file <$reference> inside inputdir <$inputdir>");
		return;
	}
	if (-e "$outputdir") {
		usage_scaffold("The outputdir <$outputdir> already exists, please choose a different outputdir");
		return;
	}
	if (!-e "$markersfile") {
		usage_scaffold("Can't find markers file <$markersfile>");
		return;
	}
	if($flow+$cytof > 1) {
		usage_scaffold("These two parameters [--flow, --cytof] can not be used jointly, please choose one of them (or none for default options)");
		return;
	}
	if($downsample ne ''){
		if($all eq '1' && $downsample ne ''){
			usage_scaffold("These two parameters [--all, --downsample NUM] can not be used jointly, please choose one of them (or none for default options)");
        	      	return;
		}
		if(!isnum($downsample) || ($downsample < 500 || $downsample > 100000)){
			usage_scaffold("<$downsample> is not a valid downsample numebr. Please insert a number between 500 and 100,000");
			return;
		}	
	}

	my $run="$ENV{CYTOFPIPE_HOME}/code/run_cytofpipe.sh "."@ARGV";
	system("$run");

}


sub parse_citrus {

	my @args=@ARGV;
        my $command=shift(@args);


	my $inputdir=''; my $outputdir=''; my $markersfile='';
	my $conditions="";
	my $flow=0;
	my $cytof=0;
	my $all='';
	my $medians='';
	my $downsample='';
	GetOptionsFromArray (
	    \@args,
	    "i=s" => \$inputdir,
	    "o=s"   => \$outputdir,
	    "m=s"   => \$markersfile,
	    "cond=s"   => \$conditions,
	    "flow"   => \$flow,
	    "cytof"   => \$cytof,
	    "all"   => \$all,
	    "downsample=i"   => \$downsample,
	    "medians=s" => \$medians,
	    "<>"   => \&print_citrus,
	) or die "\n";

	if ($inputdir eq '' || $outputdir eq '' || $markersfile eq '' || $conditions eq '') {
		usage_citrus("Please check that you are providing a inputdir (-i), outputdir (-o), markersfile (-m) and conditions file (--cond)");
		return;
	}		

	if (!-e "$inputdir") {
		usage_citrus("Can't find directory with fcs files <$inputdir>");
		return;
	}		
	if (!-e "$conditions") {
		usage_citrus("Can't find the conditions file <$conditions>");
		return;
	}
	if (-e "$outputdir") {
		usage_citrus("The outputdir <$outputdir> already exists, please choose a different outputdir");
		return;
	}
	if (!-e "$markersfile") {
		usage_citrus("Can't find markers file <$markersfile>");
		return;
	}
	if($flow+$cytof > 1) {
		usage_citrus("These two parameters [--flow, --cytof] can not be used jointly, please choose one of them (or none for default options)");
		return;
	}
	if($downsample ne ''){
		if($all eq '1' && $downsample ne ''){
			usage_citrus("These two parameters [--all, --downsample NUM] can not be used jointly, please choose one of them (or none for default options)");
        	      	return;
		}
		if(!isnum($downsample) || ($downsample < 500 || $downsample > 100000)){
			usage_citrus("<$downsample> is not a valid downsample numebr. Please insert a number between 500 and 100,000");
			return;
		}	
	}
	if($medians ne ''){
		if (!-e "$medians") {
			usage_citrus("Can't find the file <$medians> with the markers chosen to estimate medians");
			return;
		}
	}

	my $run="$ENV{CYTOFPIPE_HOME}/code/run_cytofpipe.sh "."@ARGV";
	system("$run");

}


sub print_usage {

	my @a=@_;

	my $usage0="";
	my $usage1="Program: Cytofpipe";
	my $usage2 = "Version: 1.0";
	my $usage3 = "Contact: Lucia Conde <l.conde\@ucl.ac.uk>";
	my $usage4="";
	my $usage5="Usage:   cytofpipe <command> [options]";
	my $usage6="";
	my $usage7="Commands: --clustering  Performs preprocessing of data and clustering/visualization (based on cytofkit)";
	my $usage8="          --scaffold    Constructs scaffold maps (based on Scaffold)";
	my $usage9="          --citrus      Performs citrus analysis (based on citrus)";
	my $usage10="";

	print "$usage0\n$usage1\n$usage2\n$usage3";
	print "$usage4\n$usage5\n$usage6\n$usage7\n";
	print "$usage8\n$usage9\n$usage10\n";

	die "ERROR: Invalid argument '@a'. Please select one of the available commands\n";

}


sub print_clustering {

	my @a=@_;

	my $usage0="";
	my $usage1="Program: Cytofpipe --clustering";
	my $usage2 = "Version: 1.0";
	my $usage3 = "Contact: Lucia Conde <l.conde\@ucl.ac.uk>";
	my $usage4="";
	my $usage5="Usage:   cytofpipe --clustering -i DIR -o DIR -m FILE [options]";
	my $usage6="";
	my $usage7="Required: -i DIR	Input directory with the FCS files";
	my $usage8="          -o DIR	Output directory where results will be generated";
	my $usage9="          -m FILE	File with markers that will be selected for clustering";
	my $usage10="Options: --config FILE			Configuration file to customize the analysis";
	my $usage11="         --flow | --cyto		Flow cytometry data (transformation = autoLgcl) or Cytof data (transformation = cytofAsinh) [--cytof]";
	my $usage12="         --all | --downsample NUM	Use all events in the analysis or downsample each FCS file to the specified number of events (with no replacement for sample with events < NUM) [--downsample 10000]";
	my $usage13="         --displayAll			Display all markers in output files [NULL]";
	my $usage14="";

	print "$usage0\n$usage1\n$usage2\n$usage3\n";
	print "$usage4\n$usage5\n$usage6\n$usage7\n";
	print "$usage8\n$usage9\n$usage10\n$usage11\n";
	print "$usage12\n$usage13\n$usage14\n";

	die "ERROR: Invalid argument '@a' in --clustering mode\n";
}


sub usage_clustering {
  my $error=shift;
  die qq(
Program: Cytofpipe --clustering
Version: 1.0
Contact: Lucia Conde <l.conde\@ucl.ac.uk>

Usage:   cytofpipe --clustering -i DIR -o DIR -m FILE [options]

Required: -i DIR	Input directory with the FCS files
          -o DIR	Output directory where results will be generated
          -m FILE	File with markers that will be selected for clustering
Options: --config FILE			Configuration file to customize the analysis
         --flow | --cyto		Flow cytometry data (transformation = autoLgcl) or Cytof data (transformation = cytofAsinh) [--cytof]
         --all | --downsample NUM	Use all events in the analysis or downsample each FCS file to the specified number of events (with no replacement for sample with events < NUM) [--downsample 10000]
         --displayAll			Display all markers in output files [NULL]

ERROR: $error
);
}


sub usage_clustering_config {
  my $error=shift;
  die qq(
Program: Cytofpipe --clustering
Version: 1.0
Contact: Lucia Conde <l.conde\@ucl.ac.uk>

------------------
CONFIG file format
------------------
[ clustering ]			#-- MANDATORY FIELD, IT SHOULD BE THE FIRST LINE OF THE CONFIG FILE";

GATING = yes|no			#-- MANDATORY FIELD

TRANSFORM = autoLgcl, cytofAsinh, logicle, arcsinh or none	#-- MANDATORY FIELD
MERGE = ceil, all, min or fixed				#-- MANDATORY FIELD
DOWNSAMPLE = integer between 500 and 100000			#-- MANDATORY FIELD if MERGE = fixed or ceil

#- DimRed method (tSNE) parameters:"
PERPLEXITY = 30
THETA = 0.5
MAX_ITER = 1000

#- Clustering methods:
PHENOGRAPH = yes|no
CLUSTERX = yes|no
DENSVM = yes|no
FLOWSOM = yes|no
FLOWSOM_K = number between 2 and 50				#-- MANDATORY FIELD if FLOWSOM = YES

#- Additional visualization methods:
PCA = yes|no
ISOMAP = yes|no

#- Other:
DISPLAY_ALL = yes|no
------------------

ERROR: $error
);
}


sub print_scaffold {

	my @a=@_;

	my $usage0="";
	my $usage1="Program: Cytofpipe --scaffold";
	my $usage2 = "Version: 1.0";
	my $usage3 = "Contact: Lucia Conde <l.conde\@ucl.ac.uk>";
	my $usage4="";
	my $usage5="Usage:   cytofpipe --scaffold -i DIR -o DIR -m FILE --ref FILE [options]";
	my $usage6="";
	my $usage7="Required: -i DIR		Input directory with the FCS files and the /gated subfolder with manually gated landmark populations";
	my $usage8="          --ref FILE	Name of the reference FCS file (which should be inside the Input directory)";
	my $usage9="          -o DIR		Output directory where results will be generated";
	my $usage10="          -m FILE		File with markers that will be selected for clustering";
	my $usage11="Options: --flow | --cyto		Flow cytometry data (arcsinh transformation using asinh_cofactor = 150) or Cytof data (asinh_cofactor = 5) [--cytof]";
	my $usage12="         --all | --downsample NUM	Use all events in the analysis or downsample each FCS file to the specified number of events (with no replacement for sample with events < NUM) [--downsample 10000]";
	my $usage13="";

	print "$usage0\n$usage1\n$usage2\n$usage3\n";
	print "$usage4\n$usage5\n$usage6\n$usage7\n";
	print "$usage8\n$usage9\n$usage10\n$usage11\n";
	print "$usage12\n$usage13\n";

	die "ERROR: Invalid argument '@a' in --scaffold mode\n";

}



sub usage_scaffold {
  my $error=shift;
  die qq(
Program: Cytofpipe --scaffold
Version: 1.0
Contact: Lucia Conde <l.conde\@ucl.ac.uk>

Usage:   cytofpipe --scaffold -i DIR -o DIR -m FILE --ref FILE [options]

Required: -i DIR	Input directory with the FCS files and the /gated subfolder with manually gated landmark populations
          --ref FILE	Name of the reference FCS file (which should be inside the Input directory)
          -o DIR	Output directory where results will be generated
          -m FILE	File with markers that will be selected for clustering
Options: --flow | --cyto		Flow cytometry data (arcsinh transformation using asinh_cofactor = 150) or Cytof data (asinh_cofactor = 5) [--cytof]
         --all | --downsample NUM	Use all events in the analysis or downsample each FCS file to the specified number of events (with no replacement for sample with events < NUM) [--downsample 10000]

ERROR: $error
);
}


sub print_citrus {

	my @a=@_;

	my $usage0="";
	my $usage1="Program: Cytofpipe --citrus";
	my $usage2 = "Version: 1.0";
	my $usage3 = "Contact: Lucia Conde <l.conde\@ucl.ac.uk>";
	my $usage4="";
	my $usage5="Usage:   cytofpipe --citrus -i DIR -o DIR -m FILE --cond FILE [options]";
	my $usage6="";
	my $usage7="Required: -i DIR		Input directory with the FCS files";
	my $usage8="          --cond FILE	File indicating which samples belong to each condition";
	my $usage9="          -o DIR		Output directory where results will be generated";
	my $usage10="          -m FILE		File with markers that will be selected for clustering";
	my $usage11="Options: --medians FILE		Use medians of markers in this file as statistic choice instead of cluster abundance";
	my $usage12="         --flow | --cyto		Flow cytometry data (arcsinh transformation using asinh_cofactor = 150) or Cytof data (asinh_cofactor = 5) [--cytof]";
	my $usage13="         --all | --downsample NUM	Use all events in the analysis or downsample each FCS file to the specified number of events (with no replacement for sample with events < NUM) [--downsample 10000]";
	my $usage14="";

	print "$usage0\n$usage1\n$usage2\n$usage3\n";
	print "$usage4\n$usage5\n$usage6\n$usage7\n";
	print "$usage8\n$usage9\n$usage10\n$usage11\n";
	print "$usage12\n$usage13\n$usage14\n";

	die "ERROR: Invalid argument '@a' in --citrus mode\n";

}

sub usage_citrus {
  my $error=shift;
  die qq(
Program: Cytofpipe --citrus
Version: 1.0
Contact: Lucia Conde <l.conde\@ucl.ac.uk>

Usage:   cytofpipe --citrus -i DIR -o DIR -m FILE --cond FILE [options]

Required: -i DIR	Input directory with the FCS files
          -o DIR	Output directory where results will be generated
          -m FILE	File with markers that will be selected for clustering
          --ref FILE	File indicating which samples belong to each condition
Options: --medians FILE			Use medians of markers in this file as statistic choice instead of cluster abundance
         --flow | --cyto		Flow cytometry data (transformation = autoLgcl) or Cytof data (transformation = cytofAsinh) [--cytof]
         --all | --downsample NUM	Use all events in the analysis or downsample each FCS file to the specified number of events (with no replacement for sample with events < NUM) [--downsample 10000]

ERROR: $error
);
}


sub isnum ($) {
    return 0 if $_[0] eq '';
    $_[0] ^ $_[0] ? 0 : 1
}


sub check_config_clustering {
	my $config = shift(@_);

	my $gating;my $transform;my $merge;my $downsample;my $flowsom;my $flowsom_k;
	my $perplexity; my $theta; my $max_iter;

	local $/ = undef;
	open(INF, "$config");
	my $content = <INF>;
	my @lines = split /\r\n|\n|\r/, $content;

	my $first_line=1;
	foreach my $line(@lines){
		chomp $_;
		if($first_line == 1 && $line !~ /^\s*\[\s*clustering\s*\]\s*$/i){
			usage_clustering_config("Invalid config file. Please make sure that the first line of the config file is \"[ clustering ]\"");
			return;
		}
		if($line=~/^GATING\s*\=\s*(.*)\s*$/){
			$gating=$1;
			if($gating !~/^yes$/i && $gating !~/^no$/i){
				usage_clustering_config("\"$gating\" **************** is not a valid GATING option in <$config>. Please correct the config file and choose if you want to include atomatic gating on your analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^TRANSFORM\s*\=\s*(.*)\s*$/){
			$transform=$1;
			if($transform ne "autoLgcl" && $transform ne "cytofAsinh" && $transform ne "logicle" && $transform ne "arcsinh" && $transform ne "none"){
				usage_clustering_config("Can't recognize \"$transform\" as a valid transformation method in <$config>. Please correct the config file and choose one of the available methods (\"autoLgcl\", \"cytofAsinh\", \"logicle\", \"arcsinh\", \"none\"\) or omit the config file to run with default parameters");
				return;
			}
		}
		if($line=~/^MERGE\s*\=\s*(.*)\s*$/){
			$merge=$1;
			if($merge ne "ceil" && $merge ne "all" && $merge ne "min" && $merge ne "fixed"){
				usage_clustering_config("Can't recognize \"$merge\" as a valid merge method in <$config>. Please correct the config file and choose one of the available methods (\"ceil\", \"all\", \"min\", \"fixed\"\) or omit the config file to run with default parameters");
				return;
			}
		}
		if($line=~/^PCA\s*\=\s*(.*)\s*$/){
			my $pca=$1;
			if($pca !~/^yes$/i && $pca !~/^no$/i){
				usage_clustering_config("\"$pca\" is not a valid PCA option in <$config>. Please correct the config file and choose if you want to include PCA for cluster visualization in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^ISOMAP\s*\=\s*(.*)\s*$/){
			my $isomap=$1;
			if($isomap !~/^yes$/i && $isomap !~/^no$/i){
				usage_clustering_config("\"$isomap\" is not a valid ISOMAP option in <$config>. Please correct the config file and choose if you wnat to include ISOMAP for cluster visualization in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^PHENOGRAPH\s*\=\s*(.*)\s*$/){
			my $phenograph=$1;
			if($phenograph !~/^yes$/i && $phenograph !~/^no$/i){
				usage_clustering_config("\"$phenograph\" is not a valid PHENOGRAPH option in <$config>. Please correct the config file and choose if you want to include PHENOGRAPH as clustering method in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^CLUSTERX\s*\=\s*(.*)\s*$/){
			my $clusterx=$1;
			if($clusterx !~/^yes$/i && $clusterx !~/^no$/i){
				usage_clustering_config("\"$clusterx\" is not a valid CLUSTERX option in <$config>. Please correct the config file and choose if you want to include CLUSTERX as clustering method in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^DENSVM\s*\=\s*(.*)\s*$/){
			my $densvm=$1;
			if($densvm !~/^yes$/i && $densvm !~/^no$/i){
				usage_clustering_config("\"$densvm\" is not a valid DENSVM option in <$config>. Please correct the config file and choose if you want to include DENSVM as clustering method in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^FLOWSOM\s*\=\s*(.*)\s*$/){
			$flowsom=$1;
			if($flowsom !~/^yes$/i && $flowsom !~/^no$/i){
				usage_clustering_config("\"$flowsom\" is not a valid FLOWSOM option in <$config>. Please correct the config file and choose if you want to include FLOWSOM as clustering method in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^DOWNSAMPLE\s*\=\s*(.*)\s*$/){
			$downsample=$1;
		}
		if($line=~/^FLOWSOM_K\s*\=\s*(.*)\s*$/){
			$flowsom_k=$1;
		}
		if($line=~/^PERPLEXITY\s*\=\s*(.*)\s*$/){
			$perplexity=$1;
			if($perplexity !~/^\d+$/ || $perplexity < 5 || $perplexity > 50){
				usage_clustering_config("\"$perplexity\" is not a valid PERPLEXITY option in <$config>. Please correct the config file and choose a value between 5 and 50");
				return;
			}
		}
		if($line=~/^THETA\s*\=\s*(.*)\s*$/){
			$theta=$1;
			if($theta !~/([0-9]*[.])?[0-9]+$/ || $theta < 0 || $theta > 1){
				usage_clustering_config("\"$theta\" is not a valid THETA option in <$config>. Please correct the config file and choose a value between 0 and 1");
				return;
			}
		}
		if($line=~/^MAX_ITER\s*\=\s*(.*)\s*$/){
			$max_iter=$1;
			if($max_iter !~/^\d+$/ || $max_iter < 100 || $max_iter > 5000){
				usage_clustering_config("\"$max_iter\" is not a valid MAX_ITER option in <$config>. Please correct the config file and choose a value between 100 and 5000");
				return;
			}
		}
		if($line=~/^DISPLAY_ALL\s*\=\s*(.*)\s*$/){
			my $display=$1;
			if($display !~/^yes$/i && $display !~/^no$/i){
				usage_clustering_config("\"$display\" is not a valid DISPLAY_ALL option in <$config>. Please correct the config file and choose if you want to display all the markers in the output files \(\"YES or NO\"\)");
				return;
			}
		}
		$first_line++;
	}
	close(INF);
	if(!$gating || $gating eq ''){
		usage_clustering_config("Gating parameter not found in <$config>. Please correct the config file and enter a valid Gating option \(\"YES\" or \"NO\"\) or omit the config file to run with default parameters \(GATING = NO\)");
		return;
	}
	if(!$transform || $transform eq ''){
		usage_clustering_config("Transformation parameter not found in <$config>. Please correct the config file and enter a valid transformation method \(\"autoLgcl\", \"cytofAsinh\", \"logicle\", \"arcsinh\", \"none\"\) or omit the config file to run with default parameters \(TRANSFORMATION = arcsinh\)");
		return;
	}
	if(!$merge || $merge eq ''){
		usage_clustering_config("Merge parameter not found in <$config>. Please correct the config file and enter a valid merge method \(\"ceil\", \"all\", \"min\", \"fixed\"\) or omit the config file to run with default parameters \(MERGE = fixed, DOWNSAMPLE = 10000\)");
		return;
	}
	if($merge =~ /^fixed$/i || $merge =~ /^ceil$/i){
		if(!$downsample || $downsample eq '' || $downsample !~/^\d+$/ || $downsample < 500 || $downsample > 100000){
			if(!$downsample || $downsample eq ''){
				usage_clustering_config("Downsample parameter not found in <$config>. Please correct the config file and enter a valid size between 50 and 100000");
				return;
			}elsif($downsample !~/^\d+$/ || $downsample < 500 || $downsample > 100000){
				usage_clustering_config("Can't recognize \"$downsample\" as a valid downsample number in <$config>. Please correct the config file and choose a downsample size between 500 and 100000");
				return;
			}
		}
	}
	if($flowsom =~ /^yes$/i){
		if(!$flowsom_k || $flowsom_k eq '' || $flowsom_k !~/^\d+$/ || $flowsom_k < 2 || $flowsom_k > 50){
			if(!$flowsom_k || $flowsom_k eq ''){
				usage_clustering_config("FlowSOM_k parameter not found in <$config>. Please correct the config file and enter a valid cluster number between 2 and 50");
				return;
			}elsif($flowsom_k !~/^\d+$/ || $flowsom_k < 2 || $flowsom_k > 50){
				usage_clustering_config("Can't recognize \"$flowsom_k\" as a valid number of FlowSOM clusters in <$config>. Please correct the config file and choose a cluster number between 2 and 50");
				return;
			}
		}
	}
}
