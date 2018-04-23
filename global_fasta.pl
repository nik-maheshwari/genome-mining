
system "./ggsearch36 -w 100 -d 1 -E 0.0001 -q templates.txt matLanti_full  -m 2 > glob_alignment_e-4_templates_lanpep.txt";

use Bio::SeqIO;

open $out, ">clusters_ggsearch.txt" or die;
open $in, "glob_alignment_e-4_templates_lanpep.txt" or die;
@line=<$in>;
$k=0;
for($i=0;$i<=$#line;$i++)
{
	$new="";
	chomp($line[$i]);
	if($line[$i]=~ m/>>>/gi)
	{
	$line[$i]=~ s/>//gi;
	$line[$i]=~ s/ //gi;
	$line[$i]=~ s/^[0-9]+//gi;
	@q=split("-",$line[$i]);
	$query=$q[0];
	}
	if($line[$i]=~m/>>/gi)
	{
	push @id1,$query;
	$line[$i]=~ s/>//gi;
	@h=split(" ",$line[$i]);
	$hit{$query}=$h[0];
	$line[$i+1]=~ s/ //gi;
	@e=split(":",$line[$i+1]);
	$eval{$query}=$e[-1];
	chomp($eval{$query});
	$line[$i+2]=~ s/\(//gi;
	@iden=split(" ",$line[$i+2]);
	$idperc{$query}=$iden[4];
	$sim{$query}=$iden[6];
	@al1=split(" ",$line[$i+5]);
	@al2=split(" ",$line[$i+6]);
	for($j=0;$j<length($al1[1]);$j++)
	{
		$c=substr($al1[1],$j,1);
		if(substr($al2[1],$j,1) eq ".")
		{$new=$new.$c;}
		else
		{
		$new=$new.(lc $c);
		}
	}
	$new=~ s/-//gi;
	$half=substr($new,length($new)/2);
	$chom{$query}=()=$half=~ m/\p{Lowercase}/g;
	$tempmod{$query}=$new;
	$totC{$query}=()=$new=~ m/c/gi;
	$unmatchC{$query}=()=$new=~ m/c/g;
	$k++;
	}
	
}

close $in;

open $in, "cluster.txt" or die;
@line1=<$in>;

$io=Bio::SeqIO->new(-file=>"templates.txt", -format=>"fasta");
while($prot=$io->next_seq)
{
	$id=$prot->id;
	chomp($id);
	@y=split(/\|/,$id);

	$seq1=$prot->seq;
	chomp($seq1);
	foreach $l (@line1)
	{
	chomp($l);
	@da=split("\t",$l);
	@xyz=split(/\|/,$da[-6]);
	if($da[40] eq "" && $da[-1] eq $seq1 && $xyz[-2] eq $y[-2])
	{
	$l=$l."\t".$tempmod{$id}."\t".$hit{$id}."\t".$idperc{$id}."\t".$sim{$id}."\t".$eval{$id}."\t".$totC{$id}."\t".$unmatchC{$id}."\t".$chom{$id}."\n" if $tempmod{$id} ne "";
	}
	else
	{$l=$l."\n";}
	}
}
print $out @line1;

