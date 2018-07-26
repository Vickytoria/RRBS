#!/usr/bin/perl -w
no warnings 'uninitialized';
use List::Util qw(first);


my @alphabet = ("a".."z");


#use the following variables to go to the corresponding "if"
#for example, if $usergroup=2 and $userreplic=3, we expect R,S;Ra,Rb,Rc;Sa,Sb,Sc
#and if $usergroup=3 and $userreplic=2, we expect R1,R2,R3;R1a,R1b;R2a,R2b;etc..


### SUMMARIZING, THE FOLLOWING CHUNK DECIDES WHAT SEQUENCE OF SAMPLE IDS WE SHOULD EXPECT BASED ON THE NUMBER OF GROUPS AND THE NUMBER OF REPLICATES ###
print "How many groups are there included in this analysis? (this program can only handle up to 4 groups)\n";
my $usergroups=<STDIN>;
chomp $usergroups;

my @totalsamp;
my @Namples;

if ($usergroups eq "2"){

	@Namples=("R","S");

	print "How many biological replicates in the first group?\n";
	my $Replic=<STDIN>;
	print "Your answer is $Replic";
	chomp $Replic;
	my @Ramples= ("R".$alphabet[0].."R".$alphabet[$Replic-1]);
#	print @Ramples;

	print "and how many biological replicates in the second group?\n";
	my $Seplic=<STDIN>;
	print "Your answer is $Seplic";
	chomp $Seplic;
	my @Samples= ("S".$alphabet[0].."S".$alphabet[$Seplic-1]);
#	print @Samples;
	push (@totalsamp, @Ramples);
	push (@totalsamp, @Samples);	

} else {
	
	@Namples=("R"."1".."R".$usergroups);
#	print @Namples;

	print "How many biological replicates are there in each group? (this only works if the number of replicates is the same for all groups)\n";
	my $Neplic=<STDIN>;
	print "Your answer is $Neplic";
	chomp $Neplic;
	my @Namplesneplic= ($alphabet[0]..$alphabet[$Neplic-1]);
	
	foreach my $elem (@Namples){
		foreach my $nelem (@Namplesneplic){
			push @totalsamp, "$elem$nelem";
						
		}
#	push @samples;

	}

#	print @samples;

}

#### At this point, @Namples is the complete list of group identifiers, ####
#### and @totalsamp contains all the variables we should expect for the propMeth column and in the right order ####
print "The groups are:@Namples\n";
my $n=@Namples;
my $groups = join (',', @Namples);
print "All the individual samples are:@totalsamp\n";

##### >>>>> NEW BIT
#my$lastid=$totalsamp[-1];
#print "$lastid\n";
##### >>>>> TILL HERE

my $m=@totalsamp;
my $totalsamples=join(',', @totalsamp);
#print $totalsamples;



## END OF INITIAL CHUNK ##




#### Now the actual comparison of what the file contains with what we expect begins ####
my $fname=$ARGV[0];
open (DATOS, $fname) || die "ERROR, el archivo de input no abre";
my @lines = <DATOS>;
my $header= shift @lines;
my $newline;


##adjust the header line to the new content##


if ($usergroups eq "2"){
	$header="#Chr,Start,End,Len,CpGs,Pr,Test,>Meth,Sample_count_R,Sample_count_S,$groups,$totalsamples,line\n";
	}else{
	$header="#Chr,Start,End,Len,CpGs,Pr,Test,>Meth,$groups,$groups,$totalsamples,line\n";
	}

$fname=~s/\.txt//;
my $fileOUT = "$fname.ORG.csv";
open (OUT, ">$fileOUT") || die "ERROR, el archivo de output no abre";
print OUT $header;

## a counter is introduced to add a column "line" with line numbers that will be used by the R script that runs stats later
my $count=0;

foreach my $l (@lines){
	
	$count++;
	
## to avoid the separation of the F statistic column with the comma, I change it here for the | symbol
	$l=~s/\((\d+)\,(\d+)\)/\($1\|$2\)/;
	my @line= split (/\t/, $l);
	
	## check for discrepancies between what the user asked for and what the file contains ##
	my $groups=()=$line[8]=~/\=/gi;
	if ($groups ne $usergroups){
		die "ERROR, the inserted number of groups may be wrong. Please revise.\nAnd careful!!! the number of replicates will not be double checked so you must make sure that is correct\n";
	}
		
	my @propmeth=split(";", $line[9]);
	#capture the contents of column propmeth
	my $totpropmeth = @propmeth;
	#counts how many different parts there are to this propmeth column (there is always the averages part, and then individual values per group part)
	
	pop @line;
	#remove the faulty column propmeth from each line
	#### REMEMBER REPLACING WITH THE NEWLY COMPLETED COLUMN ####
	
	# @propmeth contains as many elements as groups plus one.
	# @propmeth[0] contains the average value of propmeth per group (i.e., R=...,S=...; or R1=..., R2=..., R3=... etc...)
	# if all samples from a group are missing, that group will not appear in @propmeth[0].
	my @pm0= split(",", $propmeth[0]);
	shift @propmeth;
	# like this we have each element of the averages individually in a new array @pm0
	# and deleted that bit (averages) from the original array with all the contents
	
	my $previdx=0;	

	foreach my $pm (@pm0){
#	print "we are at $pm\n";
		$pm=~/(\w+)\=\d\.\d+/g;
		my $p=$1;
		my $idx = first { $Namples[$_] eq $p } 0..$#Namples;
#		print "The index is $idx\n";
		
		$diff=$idx-$previdx;
#		print "the difference between the index and our position in the loop is $diff\n";
		
		if ($diff ne 0){
			foreach my $i (1..$diff){
			$newline=$newline.",";
		
			}
		}
		$newline=$newline.$pm.",";
#		print "\n";
	
		$previdx=$idx;
		$previdx++;

	}
	
	### NOW TIME TO CHECK THE INDIVIDUAL SAMPLE IDS ###
	
	
	my $rest = join (',', @propmeth);

	
#	print "This is the expected sequence: @totalsamp\n";
	
#	print "THIS IS THE REST $rest\n";

	my @newpropmeth = split (',', $rest);
#	print "AND THIS IS THE NEW PROPMETH @newpropmeth\n";
	
	$previdx=0;	


	foreach my $mp (@newpropmeth){
#	print "we are at $mp\n";
		$mp=~/(\w+)\=\d\.\d+/g;
		my $p=$1;
		my $idx = first { $totalsamp[$_] eq $p } 0..$#totalsamp;
#	print "The index of indiv sample is $idx\n";
		
		$diff=$idx-$previdx;
		
		if ($diff ne 0){
			foreach my $i (1..$diff){
			$newline=$newline.",";
		
			}
		}
		
		$newline=$newline.$mp.",";
#		print "\n";
		
		$previdx=$idx;
		$previdx++;
		
		
		$newline=~s/\n,//;	
#			if ($pm0[$namp]){print "Yes, it is there\n";}else{print "Nope, not there\n";}
	}
#

###>>> for cases where the last sample ids of the line are missing, we add extra commas
###>>> until the number of columns equals the number of expected sample ids
###>>> (that is, the number of elements in @totalsamp, saved as $m)
	
	while ($previdx < $m){
###>>>	
#		print "this is the current idx $previdx and there should be $m\n";	
#		chomp $newline;
		$newline=$newline.",";
#		print "This the the $newline\n";
		$previdx++;
###>>>		
	}





###>>> because we had removed the change of line character in line 209, we replace it.
$newline=$newline.",".$count."\n";
###>>>
#chop $newline;
#print "This is the new line:\n$newline\nNEW LINE\n";
#print "this is the line as is @line\n";



#print "$newline\n";

push (@line,$newline);




$line= join(",", @line);

$line=~s/(R|S)(\d)?(\w)?=//g;


print OUT $line;
#print @line;
$newline=();




		
}




close OUT;