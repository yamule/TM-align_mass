use strict;
use warnings;

use File::Temp qw/ tempfile tempdir /;
my $tempdirr = tempdir( CLEANUP => 1 );
my $num_cpu = 8;
#print $tempdirr;

my $chain1 = "UP000000625_83333_ECOLI/AF-P36677-F1-model_v1.pdb";
my $dir2 = "UP000000625_83333_ECOLI/";
my $filelist = "UP000000625_83333_ECOLI/pdblist_b.dat";

my $tmalign_exe = "bin/TMalign_mass.exe";
my $tmalign_options = " -infmt2 9 ";

my @tempfiles;
my @temphandle;
for(my $cc = 0;$cc < $num_cpu;$cc++){
	my $tempname = $tempdirr."/temp.$cc.dat";
	push(@tempfiles,$tempname);
	my $fh;
	open($fh,">".$tempname);
	push(@temphandle,$fh);
}

my $cou = 0;
open(IN,$filelist);
while(my $ss = <IN>){
	my $fh = $temphandle[$cou%$num_cpu];
	print $fh $ss;
	$cou++;
}
close(IN);

for(my $cc = 0;$cc < $num_cpu;$cc++){
	close($temphandle[$cc]);
}
if($cou < $num_cpu){
	$num_cpu = $cou;
}
my @tmres;
my @pids;
for(my $cc = 0;$cc < $num_cpu;$cc++){
	my $outname = $tempdirr."/tempres.".$cc.".dat";
	push(@tmres,$outname);
	my $pidd = fork();
	if($pidd){
		push(@pids,$pidd);
	}else{
		system("$tmalign_exe $tmalign_options $chain1 -dir2 $dir2 ".$tempfiles[$cc]." >".$outname);
		exit(0);
	}
}
#todo process の開始時間とかで実際そのプロセスなのか判定できるならする
foreach my $pid(@pids){
	waitpid ($pid, 0);
}

for(my $cc = 0;$cc < $num_cpu;$cc++){
	open(IN,$tmres[$cc]);
	while(my $ss = <IN>){
		print $ss;
	}
	close(IN);
}
