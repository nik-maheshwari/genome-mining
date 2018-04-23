#process the hmm.pl output. Use HMMer domain cut-off < 0.01.
#Outputs 15 files, one for each domain
#Outputs a enzyme type file, containing all domains.
#Outputs all LanC or LanD sequences in parent directory.
use Bio::SeqIO;
use Cwd qw();
$path = Cwd::abs_path();
@z=split('/',$path);	#$z[-1] is the name of the species folder

my %protein=();
my $io=Bio::SeqIO->new(-file=>"translate_final.txt", -format=>"fasta");
while(my $prot=$io->next_seq)
{
	my $pro=$prot->seq;
	$pro=~ s/\*//gi;
	$protein{$prot->id}=$pro;
}


open my $outf, ">>../results_summary.txt" or die;
open my $out_type, ">enz_type.txt" or die;
open my $out_allCD, ">>../allLanCDeh.txt" or die;
open my $outbn, ">A.txt" or die;
open my $outbc, ">B.txt" or die;
open my $outc, ">C.txt" or die;
open my $outd, ">D.txt" or die;
open my $outkin, ">Kin.txt" or die;
open my $outl, ">L.txt" or die;
open my $oute, ">E.txt" or die;
open my $outv, ">V.txt" or die;
open my $outr1, ">R.txt" or die;
open my $outr2, ">Q.txt" or die;
open my $outk1, ">K.txt" or die;
open my $outk2, ">J.txt" or die;
open my $outp1, ">P.txt" or die;
open my $outp2, ">N.txt" or die;
open my $outt, ">T.txt" or die;

#open $clus, ">>../cluster.txt" or die;
my $bn=$bc=$c=$l=$k1=$k2=$kin=$p1=$p2=$t=$r1=$r2=$d=$e=$v=0;
%dupl_bn;
%dupl_bc;
%dupl_c;
%dupl_l;
%dupl_kin=%dupl_k1=%dupl_k2=%dupl_r1=%dupl_r2=%dupl_t=%dupl_p1=%dupl_p2=%dupl_d=%dupl_v=%dupl_e;
open my $in,"totalHmm.out" or die;
while(<$in>)
{
	chomp;
	if($_ =~ m/^\#/gi)
	{}
	else
	{
	my @x=split(" ");
	my @data=split('\|',$x[0]);
	my $nc=$data[-2];
	my @data1=split('_',$x[0]);
	my $id=$data1[-1];
	my $type="";
	if($x[7]<=0.01)	#HMMer domain cut off. <0.01 is considered significant
	{

	if($x[2]=~ m/LanB_N/g && !defined $dupl_bn{$x[0]})
	{
	print $outbn $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_bn{$x[0]}++;
	$bn++;
	$type="A";
	}
	elsif($x[2]=~ m/LanB_C/g && !defined $dupl_bc{$x[0]})
	{
	print $outbc $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_bc{$x[0]}++;
	$bc++;
	$type="B";
	}
	elsif($x[2]=~ m/LanC/gi && !defined $dupl_c{$x[0]})
	{
	print $outc $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_c{$x[0]}++;
	$c++;
	$type="C";
	}
	elsif($x[2]=~ m/duf/gi && !defined $dupl_d{$x[0]})
	{
	print $outd $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_d{$x[0]}++;
	$d++;
	$type="D";
	}
	elsif($x[2]=~ m/Kinase/gi && !defined $dupl_kin{$x[0]})
	{
	print $outkin $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_kin{$x[0]}++;
	$kin++;
	$type="S";
	}
	elsif($x[2]=~ m/Lyase/gi && !defined $dupl_l{$x[0]})
	{
	print $outl $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_l{$x[0]}++;
	$l++;
	$type="L";
	}
	elsif($x[2]=~ m/Bles03/gi && !defined $dupl_e{$x[0]})
	{
	print $oute $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_e{$x[0]}++;
	$e++;
	$type="E";
	}
	elsif($x[2]=~ m/VenL/gi && !defined $dupl_v{$x[0]})
	{
	print $outv $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_v{$x[0]}++;
	$v++;
	$type="V";
	}
	elsif($x[2]=~ m/NisK_1/gi && !defined $dupl_k1{$x[0]})
	{
	print $outk1 $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_k1{$x[0]}++;
	$k1++;
	$type="K";
	}
	elsif($x[2]=~ m/NisK_2/gi && !defined $dupl_k2{$x[0]})
	{
	print $outk2 $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_k2{$x[0]}++;
	$k2++;
	$type="J";
	}
	elsif($x[2]=~ m/Nis_R_1/gi && !defined $dupl_r1{$x[0]})
	{
	print $outr1 $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_r1{$x[0]}++;
	$r1++;
	$type="R";
	}
	elsif($x[2]=~ m/Nis_R_2/gi && !defined $dupl_r2{$x[0]})
	{
	print $outr2 $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_r2{$x[0]}++;
	$r2++;
	$type="Q";
	}
	elsif($x[2]=~ m/NisP/gi && !defined $dupl_p1{$x[0]})
	{
	print $outp1 $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_p1{$x[0]}++;
	$p1++;
	$type="P";
	}
	elsif($x[2]=~ m/C39/gi && !defined $dupl_p2{$x[0]})
	{
	print $outp2 $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_p2{$x[0]}++;
	$p2++;
	$type="N";
	}
	elsif($x[2]=~ m/NisT/gi && !defined $dupl_t{$x[0]})
	{
	print $outt $x[0]."_".$x[19]."_".$x[21]."\n";
	$dupl_t{$x[0]}++;
	$t++;
	$type="T";
	}
	if ($type =~ /C|D|A|B|V/)
	{
	print $out_type $nc."\t".$type."\t".$x[19]."\t".$x[21]."\t".$id."\t".$protein{$x[0]}."\n";			#$x[19] is start, $x[21] is end, $nc is accession number
	print $out_allCD ">".$z[-1]."_".$x[0]."_".$x[19]."_".$type."\n".$protein{$x[0]}."\n";
	}
	else
	{print $out_type $nc."\t".$type."\t".$x[19]."\t".$x[21]."\t".$id."\n";}						#$x[19] is start, $x[21] is end, $nc is accession number
	}
	}
}
close $out;
close $in;


print $outf $bn."-".$bc."-".$c."-".$d."-".$kin."-".$l."-".$e."-".$v."-".$k1."-".$k2."-".$r1."-".$r2."-".$p1."-".$p2."-".$t."-".$z[-1]."\n";
#print $clus "\t\t\t\tB\tC\tD\tS\tL\tP\tT\tK\tR\n";
#print $clus $z[-1]."\tTotal\t\t\t".$b."\t".$c."\t".$d."\t".$kin."\t".$l."\t".$p."\t".$t."\t".$k."\t".$r."\n";
#close ($outf,$clus);
close ($outf,$out_type,$out_allCD,$outbn,$outbc,$outc,$outd,$outkin,$outl,$outr1,$outr2,$outk1,$outk2,$outp1,$outp2,$outt,$oute,$outv);


