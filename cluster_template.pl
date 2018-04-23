#Find clusters.
#IMPORTANT! Run this only after running all other scripts in batch.pl  
use Bio::SeqIO;
use List::Util qw[min max];
system "cp enz_type.txt enz_type_template.txt";
open my $out, ">>enz_type_template.txt" or die;
our %duplicates=%duplicates_new=%template=%idtemp=();
my $io1=Bio::SeqIO->new(-file=>"template.txt",'-format'=>'Fasta');
while($seq_dna=$io1->next_seq)

{	
	my @simple=();
	my $print_id=$print_new="";
	my $flag=$countEnz=0;
	my $pep=$seq_dna->seq;
	chomp($pep);
	my $countC=$countS=$countT=$posC=0;
	$countC=()=$pep=~ m/C/g;
	$countS=()=$pep=~ m/S/g;
	$countT=()=$pep=~ m/T/g;
	$posC=rindex($pep,'C');	#last position of C . && length($pep)>28
	if($countC>0 && (($countS+$countT)>$countC) && (($posC+1)>2*(length($pep)/3))  )	#>1 cysteine, ser+thr>cys, last cys in last third of the peptide, removes duplicate peptides && !defined $duplicates{$pep}
	{
	my $id1=$seq_dna->id;
	my @data=split('\|',$id1);
	my $new_id=$data[-2];
	my @data1=split("_",$id1);
	my $id2=$data1[-1];
	my $desc=$seq_dna->desc;
	$desc=~s/ //gi;
	my @des=split('#',$desc);
	my $temp_start=$des[1];
	my $temp_end=$des[2];
	$duplicates{$pep}++;
	$template{"$new_id"}{"$temp_start"}=$pep;
	$idtemp{"$new_id"}{"$temp_start"}=$id1;
	print $out $new_id."\tZ\t".$temp_start."\t".$temp_end."\t".$id2."\n";
	}
}
close $out;



my $i=0;

#open $in, "test.txt" or die;
#open $out, ">test.txt" or die;
open $out, ">>../cluster.txt" or die;
open $outtemp, ">>../templates.txt" or die;

use Cwd qw();
$path = Cwd::abs_path();
our @z=split('/',$path);	#$z[-1] is the name of the species folder
my %cycle=();
open my $in, "enz_type_template.txt" or die;
while(<$in>)
{
	chomp;
	if(!defined $duplicates_new{$_})
	{
	push our @data5, $_;
	$duplicates_new{$_}++;
	}
}
close $in;
# Sort by accession number and then prot start position

my @sortedName = sort { (split ' ', $a)[0] cmp (split ' ', $b)[0] || (split ' ', $a)[2] <=> (split ' ', $b)[2] }@data5;
foreach my $ele (@sortedName)
{	# NC	Type    Start  End	ID    LanC
	($q[$i],$w[$i],$e[$i],$r[$i],$id[$i],$cycle[$i])=split(" ",$ele);
	$i++;
	
}

open my $infile, "../results_summary.txt" or die;
while(<$infile>)
{
	chomp;
	if($_=~m/$z[-1]/gi)
	{
	my @y=split("-",$_);
	#@x=split("\/",$y[-1]);
	pop(@y);
	our $total= join "\t",@y;
	our $total=$total."\t".$z[-1];
	}

}
close $infile;



sub cluster_type()
{
	my $tem=$_[0];
	my $tempcount=$_[1];
	my %countnew=%{$_[2]};
	my $cl=my $cl1=$_[3];
	my $st=$_[4];
	my $en=$_[5];
	my $f=$_[6];
	my @printC=@{$_[7]};
	my $strin="";
	my ($score,$sum,$score3)=0;
	$cl1=~ s/\.//gi;
	print  $out $st."\t".$en."\t".($en-$st+1)."\t";
	@files=qw(A B C D S L E V K J R Q P N T);
	foreach $k1 (@files)
	{
	if ($countnew{$k1}>0)
	{print  $out $countnew{$k1}."\t";
	}
	else
	{print  $out "0\t";
	}
	$score++ if $countnew{$k1}>0;		#UDC
	}
	#print $tem."\n";
	@d=split("\t",$total);
	
	$i1=0;
	foreach $k1 (@files)
	{
	if($countnew{$k1}>0)
	{
	if($d[$i1]-$countnew{$k1}==0)
	{$score3++;}
	else
	{$score3+=(1/($d[$i1]-$countnew{$k1}));}	#CGR
	}
	$i1++;
	}
	@data3=split(":",$tem);
	print  $out "$total:".$idtemp{"$data3[0]"}{"$data3[1]"}."\t$cl\t$cl1\t$score\t$score3\t";
	print $out "\n" if $f==0;
	
	print $out $template{"$data3[0]"}{"$data3[1]"}."\n" if ($f>0);
	foreach $cc (@printC)
	{$strin.="&".$cc;}
	#print $idtemp{$template{"$data3[0]"}{"$data3[1]"}}."\n";
	print $outtemp ">".$z[-1].":".$idtemp{"$data3[0]"}{"$data3[1]"}."_".$cl."_".$f."\n".$template{"$data3[0]"}{"$data3[1]"}."\n" if ($f>0);
}


#Remove clusters with length=1 or end-start>50000 or identical type (say,TTT or RR. Remove them)
sub print_cluster()
{
	$flag=$u=$v=0;
	%count=();
	$clus=$_[0];
	$clusnew=uc $_[0];
	$start1=$_[1];
	$end1=$_[2];
	$len=length($clus);
	@templ=@{$_[3]};
	@protC=@{$_[4]};
	$mini=$_[5];
	$more=$_[6];

	$clusnew=~ s/Z//gi;
	$clusnew=~ s/\.//gi;
	#print $clusnew."\n";
	#print $clus."\t".$clusnew."\t";
	for(split //,$clusnew)
	{
	++$count{$_};
	$flag=1 if($count{$_}==length($clusnew));	#Removes cluster of identical type (say ZTTTZ) 
	}
	#print $flag."\n";
	if((length($clusnew)>1 && ($end1-$start1)<=50000 && $flag==0 && $count{C}>0) || $more==5)
	#if(length($clus)>1 && ($end1-$start1)<=20000 && $flag==0 && (($count{B}>0)||($count{C}>0)||($count{D}>0)||($count{S}>0)||($count{L}>0)))
	{
	#foreach $lan (@protC)
	#{
	#print $outlanC ">".$clus."_".$z[-1]."_".$start1."_".$end1."_".$u++."\n".$lan."\n"; #if($count{B}>0||$count{A}>0||$count{D}>0||$count{L}>0||$count{E}>0||$count{V}>0||$count{S}>0);	#Print LanC for different dehydration clusters
	#print $mintemp $clus."_".$z[-1]."_".$start1."_".$end1."_".$v++."\t".$mini."\n";
	#print $outlanCB ">".$clus."_".$z[-1]."_".$start1."\n".$lan."\n" if($count{B}>0||$count{b}>0);
	#print $outlanCD ">".$clus."_".$z[-1]."_".$start1."\n".$lan."\n" if($count{D}>0);
	#print $outlanCL ">".$clus."_".$z[-1]."_".$start1."\n".$lan."\n" if($count{L}>0);
	#print $outlanCE ">".$clus."_".$z[-1]."_".$start1."\n".$lan."\n" if($count{E}>0);
	#print $outlanCV ">".$clus."_".$z[-1]."_".$start1."\n".$lan."\n" if($count{V}>0);
	#}

	if(@templ)
	{
	
	$y=1;
	foreach $t (@templ)
	{
	&cluster_type($t,$y,\%count,$clus,$start1,$end1,$y,\@protC);
	$y++;
	}
	}
	else
	{
	$t="";
	&cluster_type($t,$y,\%count,$clus,$start1,$end1,0,\@protC);
	}
	}

}

@store=@lanc=();
$lost=0;
$min=10000;
$start=$e[0];
$nc=$q[0];
$end=$r[0];
$cluster=$w[0];
$pod=$r[0] if $w[0] =~ m/B|C|A|D|L|E|V/;
$idc=$id[0];
$h=$z=$l=0;
if($cluster eq "Z")
{$store[0]=$nc.":".$start;$h=1;$z++;$temp=$e[$j];}
if($cluster eq "C")
{$lanc[0]=$cycle[0];$l=1}


for($j=1;$j<$i;$j++)
{
	
	if($q[$j] eq $nc)
	{
		#if($end-$start > 50000 && !($w[$j] =~ m/B|C|A|D|L|E|V/))
		#{&print_cluster($cluster,$start,$end,\@store,\@lanc,$min,5) if(length($cluster)>1);}
		if(($e[$j]-$end) < 10000)
		{	
			
			$point="";
			$gap=0;
			$gap=$id[$j]-$idc;	#number of non-lan proteins between the two lan proteins
			for($pt=1;$pt<$gap;$pt++)
			{
			$point.=".";
			}
			$pod=$r[$j] if $w[$j] =~ m/B|C|A|D|L|E|V/;
			if($e[$j]-$end<0)
			{$w[$j]=lc $w[$j];}
			$cluster="$cluster$point$w[$j]";
			if($w[$j] eq "Z" || $w[$j] eq "z")
			{$store[$h]=$q[$j].":".$e[$j];$h++;$temp=$r[$j];}
			if($w[$j] eq "C" || $w[$j] eq "c")
			{$lanc[$l]=$cycle[$j];$l++;}
			$idc=$id[$j];
			#print $j."\t".$cluster."\t".$end."\n";
			#$start=$e[$j];
			$end=$r[$j];
			$nc=$q[$j];
			if($w[$j] =~ m/B|C|A|D|L|E|V/i)
			{$min=min(abs($e[$j]-$temp),$min);}
			elsif($w[$j]=~ m/Z/i)
			{$min=min(abs($e[$j]-$pod),$min);}
			else {}
			$lost=5 if $cluster=~ m/B|C|A|D|L|E|V/i;
		}
		else
		{
		
		&print_cluster($cluster,$start,$end,\@store,\@lanc,$min,$lost) if(length($cluster)>1);
		$cluster=$w[$j];
		$start=$e[$j];
		$end=$r[$j];
		$nc=$q[$j];
		$lost=0;
		$min=10000;
		$pod=$r[$j] if $w[$j] =~ m/B|C|A|D|L|E|V/;
		$idc=$id[$j];
		@store=@lanc=();
		$h=$z=$l=0;
		if($w[$j] eq "C")
		{$lanc[$l]=$cycle[$j];$l++;}
		if($cluster eq "Z")
		{$store[0]=$nc.":".$start;$h=1;$z++;$temp=$e[$j];}
		}
	}
	else
	{
	&print_cluster($cluster,$start,$end,\@store,\@lanc,$min,$lost) if(length($cluster)>1);	
	$nc=$q[$j];
	$start=$e[$j];
	$end=$r[$j];
	$cluster=$w[$j];
	$idc=$id[$j];
	$lost=0;
	$min=10000;
	$pod=$r[$j] if $w[$j] =~ m/B|C|A|D|L|E|V/;
	@store=@lanc=();
	$h=$z=$l=0;
	if($w[$j] eq "C")
	{$lanc[$l]=$cycle[$j];$l++;$cluster=$cluster."C";}
	if($cluster eq "Z")
	{$store[0]=$nc.":".$start;$cluster="Z";$h=1;$z++;$temp=$e[$j];}
	}
}
	#print @lanc;	
	&print_cluster($cluster,$start,$end,\@store,\@lanc,$min,$lost) if(length($cluster)>1);	#print last cluster
	#print $out "\n";






