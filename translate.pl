#Translate genomic sequence to protein sequences, using Prodigal
use Bio::SeqIO;
opendir (DIR, '.') or die $!;

while(my $file = readdir(DIR))
{
	if($file =~ m/$\.fasta/)
	{
	my $seqlen=0;
	my $temp="";
	my $io1=Bio::SeqIO->new(-file=>$file, -format=>'Fasta');
	while(my $p=$io1->next_seq)
	{
	my $seq1=$p->seq;
	$seqlen+=length($seq1);
	}
	if($seqlen<100000)
	{$temp="-p meta";}
	my $cmd = join('', '../prodigal -q ',$temp,' -i ', $file, ' -a translate_prod_', $file, ' >> log.txt');		#running in single mode or meta mode (if overall genomic sequence is < 100000bp)
	system ($cmd);
	print $file." translated\n";
	}
	
}
system('cat translate_prod_* > translate_final.txt');
open my $out_tem, ">template.txt" or die;
open my $out_enz, ">enzyme.txt" or die;
my $io=Bio::SeqIO->new(-file=>'translate_final.txt', -format=>'Fasta');
while(my $seq_dna=$io->next_seq)
{
	my $seq=$seq_dna->seq;
	my $id=$seq_dna->id;
	my $desc=$seq_dna->desc;
	chomp($seq);
	$seq=~ s/\*//gi;
	my $len=length($seq);
	if($len<=100)		#Precursor length <=100 aa.
	{
	print $out_tem ">$id $desc\n$seq\n";
	}
	else			#Other sequences are searched against Pfam HMM database.
	{
	print $out_enz ">$id $desc\n$seq\n";
	}
}
close ($out_tem,$out_enz);
