#to allow for multiple fasta file for the same strain, have strain_file1.fasta, strain_file2.fasta, etc. 
#No _ in strain name
#have no important folder or .txt file in the location
#run from the diretory where all .fasta files are present
use strict;
use warnings;
use Cwd;
use File::chdir;

system 'rm -rf ./*/';
system 'rm -rf ./*.txt';
opendir (DIR, '.') or die $!;
while(my $file = readdir(DIR))
{
	if($file =~ m/\.fasta$/)
	{
	my @x=split("_",$file);			
	if(-d $x[0])
	{}
	else
	{
	my $cmd = join('', 'mkdir ', $x[0]);
	system ($cmd);
	}
	my $cmd1 = join('', 'cp ', $file , ' ', $x[0]);
	system($cmd1);
	}
}
closedir(DIR);
system 'for d in */; do cp batch.pl translate.pl hmm.pl process_hmm.pl cluster_template.pl "$d"; done;';	#copies .pl scripts to each folder in the current directory


#execute batch.pl inside each species folder

while (my $node = glob '*' )
{
    	next unless -d $node;
    	local $CWD = $node;
    	chdir($node);
	system("perl batch.pl");
}

#calculates global identity of predicted precursors to known lanthipeptides

system "perl global_fasta.pl";

#adds predicted bridge to each lanthipeptide. calc min distance to cyclase and Dehydratase

system "perl addBridge_lanCDeh.pl";

#clusters created. file: clusters_ggsearch_bridge_lanCDeh.txt


