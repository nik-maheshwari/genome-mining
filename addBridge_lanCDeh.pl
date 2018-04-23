#Add predicted bridge to each template
#For each template in cluster, calculates min distance to cyclase and Dehydratase. If no Deh in cluster, looks outside the cluster.


use List::Util qw[min max];
use strict;
use warnings;
use lib "./";		#change to location where Bridge.pm is stored. Default is current directory.
use Bridge qw(&brid);
my $parent = "/home/nik/Documents/code";	#change to your current directory. Give full path.
open my $out, ">clusters_ggsearch_bridge_lanCDeh.txt" or die;

print $out "Start	End	Length	A	B	C	D	S	L	E	V	K	J	R	Q	P	N	T	Tot_A	Tot_B	Tot_C	Tot_D	Tot_S	Tot_L	Tot_E	Tot_V	Tot_K	Tot_J	Tot_R	Tot_Q	Tot_P	Tot_N	Tot_T	Species	Cluster_gaps	Cluster	UDC	CGR	Template	Mod_Temp	Lanthipeptide	iden	sim	E-val	TotalC	UnmatchC	Unmatch2ndHalf	Bridge	AccessionID	PepID	PepLen	MinC	MinLanC	MinD	MinLanD\n";

open my $in, "clusters_ggsearch.txt" or die;
while(<$in>)
{

	my @lanc=my @lanD=my @x=my @temp=my @y=my @w=my @w1=my @z=();
	my %close=my %closeD=();
	my $i=my $st=my $en=my $len=0;
	my $min=my $minD=10000000;
	
	my $line=$_;
	chomp($line);
	@x=split("\t",$line);
	my $seq=$x[38];		#predicted lanthipeptide
	my ($bridge,$pos,$temp1)= &brid($seq);	#get topology. Only $bridge required.
	my $countC=()=$seq=~ m/C/g;	#count Cysteine

	my $org1=$x[33];
	@z=split(":",$org1);
	my $org=$z[0];
	@temp=split(/\|/,$z[1]);
	my $protID=$temp[-1];
	$protID=~ s/_//g;
	my $nc=$temp[-2];

	my $path = $parent."/".$org;
	print $path;
	chdir($path) or die;

	open my $in1, "enz_type_template.txt" or die;
	while(<$in1>)
	{
	chomp;
	@y=split("\t",$_);
	if($y[4] eq $protID)
	{
	$st=$y[2];
	$en=$y[3];
	}
	if($y[1] =~ /C/ and $y[0] eq $nc)		#collecting all cyclase sequence location
	{
	push @lanc, $_;
	}
	elsif($y[1] =~ /D|A|B|V/ and $y[0] eq $nc)		#collecting all Dehydratase sequence location
	{
	push @lanD, $_;
	}
	else{}
	}
		
	close $in1;
	#calculating closest LanC
	foreach my $f (@lanc)
	{
	@w=split("\t",$f);
	$min=min(abs($st-$w[3]),abs($en-$w[2]),$min);
	$close{$min}=$w[5] if (abs($st-$w[3])<=$min || abs($en-$w[2])<=$min);
	}
	#calculating closest dehydratase
	foreach my $g (@lanD)
	{
	@w1=split("\t",$g);
	$minD=min(abs($st-$w1[3]),abs($en-$w1[2]),$minD);
	$closeD{$minD}=$w1[5] if (abs($st-$w1[3])<=$minD || abs($en-$w1[2])<=$minD);
	}
	$len=($en-$st+1)/3 if $st!=0;
	if($x[38] eq "")		#no predicted lanthipeptide in the cluster
	{$min=$minD="NA";}
	if(!$x[39])		#no homology to any known lanthipeptide
	{
	print $out $line."\tNA\tNA\t0\t0\t1\t$countC\t$countC\t".length($seq)."\t".$bridge."\t".$nc."\t".$protID."\t$len\t$min\t".$close{$min}."\t".$minD."\t".$closeD{$minD}."\n";
	}
	else
	{
	print $out $line."\t".$bridge."\t".$nc."\t".$protID."\t$len\t$min\t".$close{$min}."\t".$minD."\t".$closeD{$minD}."\n";
	}

}
close $in;
close $out;
