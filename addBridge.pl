use strict;
use warnings;
use Bio::SeqIO;
use lib "./"; #change this path to location of Bridge.pm
use Bridge qw(&brid);

open my $out, ">peptide_bridge.txt" or die;
print $out "ID\tSequence\tBridge\n";
my $io=Bio::SeqIO->new(-file=>"peptide.fasta",'-format'=>'Fasta'); #Peptides in FASTA format
while(my $pep=$io->next_seq)
{
	my $id=$pep->id;
	my $seq=$pep->seq;
	my ($bridge,@pos,%temp1)= &brid($seq);
	print $out $id."\t".$seq."\t".$bridge."\n";
}

close $out;
