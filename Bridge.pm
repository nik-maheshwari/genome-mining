#perl module to calculate the bridging topology of a lanthipeptide precursor.
#returns topology code (ex SCSC) and bridging position (24,28,30,34)
#to use
#use lib "/home/nikunj/module/project/"; #change this path to location of Bridge.pm
#use Bridge qw(&brid);
#($bridge,@pos,%temp1)= &brid($seq);

package Bridge;

#use strict;
#use warnings;
no warnings 'experimental::smartmatch';
use Exporter qw(import);
our @EXPORT_OK = qw(brid);
use Bio::SeqIO;
use List::MoreUtils qw(firstidx);
	
sub brid()
    {
	$protein=$cton=$temp="";
	$len=$firstC=0;
	@lab=@positions=@labvar=();
	%con=%taken=%label=();
	$protein=shift;
	$len=length($protein);
    	$cton=reverse $protein;
	#----------get the location of C in labionin in reverse pep, and search if that C matches with CS. If yes, rename CS to CL and then to LLC----#
	while($cton=~ m/C[^C]{2,5}S[^C]{2}S/g)	#Real Lab S[^C]{2}S[^C]{2,5}C
	{
	push @lab, ($-[0]+1);
	}
	#foreach $g (@lab)
	#{print $g."\t";}

	#-------------------------------#
    	$firstC=index($cton,'C');

    	while($firstC!=-1)
    	{
	$flag=0;
    	$sub=substr($cton,$firstC);
    	$countS=()=$sub=~ m/S/g;
	$countT=()=$sub=~ m/T/g;
	
    	if(($countS+$countT)>0)
    	{
   		#Predicting lanthionine and labionin bridges (C-S/C-T, min 3 aa gap)
    		$S=index($cton,"S",$firstC+3);
		$T=index($cton,"T",$firstC+3);
		if(($firstC+1) ~~ @lab)
		{
		$con{$S+1}="S";
		$taken{$S}++;
		$SS=index($cton,"S",$S+3);
		$con{$SS+1}="S";
		$label{$S+1}=$firstC+1;
		$label{$SS+1}=$S+1;
		$taken{$SS}++;
		$flag=1;
		}
		$S=index($cton,"S",$S+1) if defined $taken{$S};
    		$T=index($cton,"T",$T+1) if defined $taken{$T};
	if($S==-1 && $T==-1)
	{$offset="";}
    	elsif($S==-1)
    	{$offset=$T;}
    	elsif($T==-1)
    	{$offset=$S;}
    	else
    	{
        if($S<$T) #S is likely to be bridged with C.
        {
        $offset=$S; #position of S on the reverse sequence
        }
        else
        {
        $offset=$T; #position of T on the reverse sequence
        }
    	}
    	
    	}
	$taken{$offset}++;
	$con{$firstC+1}="C";
	if ($offset ne "" and $flag==0)
	{$con{$offset+1}="S" ;
	$label{$offset+1}=$firstC+1;}
	$firstC=index($cton,'C',$firstC+1);
   	$offset="";
   	}
	foreach $f (sort {$a <=> $b} keys %con)
	{
	if($f ~~ @lab)
	{$con{$f}="L";}
	$temp.=$con{$f};
	push @positions, $len-$f+1;
	}
	$bridge=reverse $temp;
	@positi=reverse @positions;
	$bridge=~ s/SSL|SL/LLC/g;		#LLC denotes Lab
	return ($bridge,\@positi,\%label);
    }
1;


