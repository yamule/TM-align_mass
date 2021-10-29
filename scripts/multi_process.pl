use strict;
use warnings;

use File::Temp qw/ tempfile tempdir /;
use Getopt::Long qw(GetOptions);

use File::Basename;
my $dirname = dirname(__FILE__);

my $tempdirr = tempdir( CLEANUP => 1 );
my $num_cpu = 8;
my $tmscore_threshold = 0.4;
#print $tempdirr;

my $chain1;
my $dir2;
my $filelist;

my $helpmessage = "Usage: $0 --query query_pdb_file --dir root_dir_for_subject_files --list subject_file_name_list [--num_threads number_of_threads]  [--tmscore_threshold threshold_for_tmscore(normalized by length of Chain_1)] \n";
GetOptions('query=s' => \$chain1,'dir=s' => \$dir2, 'list=s' => \$filelist, 'num_threads:i' => \$num_cpu, '--tmscore_threshold:f' => \$tmscore_threshold) or die $helpmessage;


my $tmalign_exe = $dirname."/../bin/TMalign_mass.exe ";
my $tmalign_options = " -outfmt 2 ";

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
open(IN,$filelist) or die "Can not open $filelist.";
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

my @res;
my $head = "";
for(my $cc = 0;$cc < $num_cpu;$cc++){
	open(IN,$tmres[$cc]);
	my $head_ = <IN>;
	if($cc == 0){
		$head = $head_;
	}
	while(my $ss = <IN>){
		if($ss =~ /^Total CPU time/){
			next;
		}
		my @ptt = split(/\t/,$ss);
		if($#ptt < 2){
			next;
		}
		if($ptt[2] > $tmscore_threshold){
			my %tmp;
			$tmp{"tmscore"} = $ptt[2];
			$tmp{"line"} = $ss;
			push(@res,\%tmp);
		}
		#print $ss;
	}
	close(IN);
}

my @sortedd = sort{${$b}{"tmscore"} <=> ${$a}{"tmscore"}}@res;
print $head;
foreach my $ss(@sortedd){
	print ${$ss}{"line"};
}