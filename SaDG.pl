#!/usr/bin/perl

# To-Do list::
# improvement #1 : create more rules for GTF format.

use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $version="1.1";
my $tmp_DIR="/tmp";
#my $Nthreads=16;
# Program paths defined by the user:
my $HTseqcount_PATH="/local/jperales/Soft/NGS/HTSeq-0.5.4p5/scripts/htseq-count";
my $samtools_PATH="/local/jperales/Soft/NGS/samtools-1.2/samtools";
#


#########
#---- FUNCTIONS
# usage
sub usage($){
	my $version=$_[0];
	print STDERR "\n----\n";
	print STDERR "SaDGcurve :: Saturation Curve for Detected Genes from RNAseq data","\n";
	print STDERR "Usage ::","\n";
	print STDERR "perl SaDG.pl -G <genes.gtf> \\ \n",
		"-o <output_prefix>"," \\ \n",
		"-RD <Million Reads>"," \\ \n",
		"-by <Million Reads>"," \\ \n",
		"-R <repeats of by-step>"," \\ \n",
		"-S sampleName1:<input1.bam,input2.bam,inputN.bam>"," \\ \n",
		"-S sampleName2:<input1.bam,input2.bam,inputN.bam>"," \\ \n",
		"-S sampleNameN:<input1.bam,input2.bam,inputN.bam>","\n";
	exit 1;
}

# parameters
sub get_params($$){
	my @args=@{$_[0]};
	my $v=$_[1];

	my ($GTF,$output,$RD,$BY,$REP);
	my %samples;
	for(my $i=0;$i<scalar(@args);$i++) {
		if($args[$i] eq "-G") {
			$i++; 
			$GTF=$args[$i];
			check_GTF($GTF);
			next;
		} elsif($args[$i] eq "-o") {
			$i++;
			$output=$args[$i];
			# Future: improv #2 : check for dirname
			next;
		} elsif($args[$i] eq "-RD") {
			$i++;
			$RD=$args[$i];
			if($RD !~ /[0-9]{1,2}/){
				print STDERR "ERROR : RD is not a correct number ",
				"(1 or 2 digits).","\n";
				exit 1;
			}
			next;
		} elsif($args[$i] eq "-by") {
			$i++;
			$BY=$args[$i];
			if($BY !~ /[0-9]/) {
				print STDERR "ERROR : by must be a number","\n";
				exit 1;
			} elsif($BY > $RD) {
				print STDERR "ERROR : by must be less than Read depth","\n";
				exit 1;
			}
			next;
		} elsif($args[$i] eq "-R") {
			$i++;
			$REP=$args[$i];
			if($REP !~ /^[0-9]+$/){
				print STDERR "ERROR : repeats must be a integer number","\n";
				exit 1;
			}
			next;
		} elsif($args[$i] eq "-S") {
			$i++;
			if($args[$i] =~ /^([a-zA-Z0-9_]+):(.*)/) {
				my $sname=$1;
				my $BAMs=$2;
				if(exists($samples{$sname})){
					print STDERR "ERROR: the same Sample name",
					" was used twice:'$sname'","\n";
					exit 1;
				}
				foreach my $bam(split(",",$BAMs)){
					check_BAM($bam);
				}
				$BAMs=~ s/,/ /g;
				$samples{$sname}=$BAMs;
			}
			next;	
		} else {
			print STDERR "ERROR : Parameter not recognized:",
			"'",$args[$i],"'\n";
			usage($v);
		}
	}
	return($GTF,$output,$RD,$BY,$REP,\%samples);
}

sub check_GTF($){
	my $gtf=$_[0];
	if(-e $gtf){
		open(IN,"$gtf") or die("Cannot open $gtf : $!");
		while(my $line=<IN>){
			next if($line =~ /^#/);
			my @fds=split("\t",$line);
			if(scalar(@fds)!=9) {
				print STDERR "ERROR: BAD GTF format at '$gtf'","\n";
				exit 1;
			}
			# Future: improvement #1
		}
		close(IN);

	} else {
		print STDERR "ERROR : input file does not exist. ";
		print STDERR "'$gtf' does not exist.","\n";
		exit 1;
	}
}

sub check_BAM($){
	my $bam=$_[0];
	unless(-e $bam){
		print STDERR "ERROR : input file does not exist:";
		print STDERR "'$bam'","\n";
		exit 1;
	}
}
#---


#--- pipeline
my ($GTF,$output,$RD,$BY,$REP,$ref_samples) = get_params(\@ARGV,$version);

my %samples=%{$ref_samples};

# Get outputs
my $output_tsv=$output."_table.tsv";
my $logerr_txt=$tmp_DIR."/SaDG.logerr";

# Transform to million scale
$RD=$RD*1000000;
$BY=$BY*1000000;
my @output_table;
foreach my $sample(keys(%samples)){
	my $BAMs=$samples{$sample};
	my $shufBAM_tmp=$tmp_DIR."/SaDG_".$sample;

	my $Sys_cmd=$samtools_PATH." merge - ".$BAMs." | ".
		$samtools_PATH." bamshuf - ".
		$shufBAM_tmp.
		" 2> ".$logerr_txt;
system($Sys_cmd);
	$shufBAM_tmp.=".bam";

	$Sys_cmd=$samtools_PATH." flagstat ".$shufBAM_tmp.
		" 2>> ".$logerr_txt;
	
	my @tmp=`$Sys_cmd` or die("ERROR: flagstat couldnt :$! \n");
	my $tmp_ln1= $tmp[0];
	my ($Nreads_pass,$Nreads_failed) = ($tmp[0] =~ /^([0-9]+) \+ ([0-9]+)/);
	my $Nreads = $Nreads_pass+$Nreads_failed;

	# Fake line
#	push(@output_table,$sample."\t0\t0");
	for my $iter (1 .. $REP) {
		print STDOUT $sample,"\t","0","\t","0","\t","0","\t","0","\n";
	}

	for(my $subRD=$BY;$subRD <= $RD;$subRD+=$BY) {
		for my $iter (1 .. $REP){
		my $last_loop=0;
		# Check if subRD is greater than total reads
		if($subRD >  $Nreads) {
			$last_loop=1;
			$subRD=$Nreads;
		}

		# Generate the downsampling ids
		my @arrIDs=shuffle(1 .. $Nreads);
		@arrIDs=splice(@arrIDs,0,$subRD);
		# Sorting by descending numeric order, bc pop is more efficient than unshift
		@arrIDs= sort {$b <=> $a} @arrIDs;

		# Downsampling file
		my $down_sam=$tmp_DIR."/SaDG_".$sample."_".$subRD."_".$iter.".sam";
		open(OUT,">$down_sam") or die("CAnnot write $down_sam : $!");

		# Read the mother
		open(SAM,"$samtools_PATH view $shufBAM_tmp |") or 
			die("Cannot open $shufBAM_tmp: $!");
		my $reads_grab=0;
		while(scalar(@arrIDs)!=0){
			my $new_ID=pop(@arrIDs);
			while(my $readAln=<SAM>) {
				$reads_grab++;
				if($reads_grab==$new_ID){
					print OUT $readAln;
					last;
				}
			}	

		}
		close(SAM);
		# Gene Expression Quantification
		$Sys_cmd = "python ".$HTseqcount_PATH." --stranded "."no".
			" --mode "."intersection-nonempty".
			" --minaqual "."0".
			" --type "."exon".
			" --idattr "."gene_id".
#			" - ".$GTF." 2> ".$tmp_DIR."/SaDG.logerr";
			" ".$down_sam." ".$GTF." 2>> ".$tmp_DIR."/SaDG.logerr";
			
			my @lines=`$Sys_cmd` or 
			die("Could not exec foo: $! \n Last CMD '$Sys_cmd' \n");

			my $lib_size=0;
			foreach my $line(@lines) {
				chomp($line);
				my @fds=split("\t",$line);
				next if(scalar(@fds)>2);
				next if($fds[0] =~ /^__/); # No feature
				$lib_size=$lib_size+$fds[1];
			}
			my $counts_cnt=0;
			my $counts10_cnt=0;
			my $cpm_cnt=0;
			foreach my $line(@lines){
				chomp($line);
				my @fds=split("\t",$line);
				next if(scalar(@fds)>2);
				next if($fds[0] =~ /^__/); # No feature

				# Add+1 if 10 reads for a gene
				$counts10_cnt++ if($fds[1] >= 10);

				# Add+1 if 1 reads for a gene
				$counts_cnt++ if($fds[1] >= 1);

				# Add+1 if 1 cpm for a gene
				my $cpm=1000000*$fds[1]/$lib_size;
				$cpm_cnt++ if($cpm >= 1);				
				#print $line,"\n" if($fds[1] >= 1);
			}
#		push(@output_table,$sample."\t".$subRD."\t".$cnt);
		print STDOUT $sample,"\t",$subRD,
				"\t",$cpm_cnt,
				"\t",$counts10_cnt,
				"\t",$counts_cnt,"\n";
	
		# Remove the temporary bam file for this iteration
		unlink($down_sam);

		# Last loop?
		last if($last_loop==1);
		}
	}
	unlink($shufBAM_tmp);
}

exit;
