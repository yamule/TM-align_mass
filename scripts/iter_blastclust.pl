use strict;
use warnings;


#perl iter_blastclust.pl -in samples/small_fas.fas -out clustered.dat -outdir clusout -ident 30 -cov_long 0.5 -cov_short 0.5
my $blastppath = "D:/dummy/programs/blast/ncbi-blast-2.10.1+/bin/blastp.exe ";

#my $makeblastdbpath = "D:/dummy/programs/blast/ncbi-blast-2.10.1+/bin/makeblastdb.exe"; #Ç®ÇªÇÁÇ≠ÉoÉOÇ≈ DB çÏê¨Ç™Ç≈Ç´Ç»Ç¢
my $makeblastdbpath = "D:/dummy/work/aptamer/exapps/ncbi-blast-2.6.0+/bin/makeblastdb.exe ";



sub get_seq{
	my $filename = $_[0];
	my $buff = "";
	my $name = "";
	my $desc = "";
	my @ret;
	open(IN,$filename);
	while(my $ss = <IN>){
		if($ss =~ /^[\s]*>[\s]*([^\s]+)/){
			my $tname = $1;
			my $tdesc = "";
			if($ss =~ /^[\s]*>[\s]*[^\s]+[\s]+([^\s].+)/){
				$tdesc = $1;
			}
			if(length($buff) > 0){
				my %tmp;
				$buff =~ s/[\s]//g;
				$tmp{"seq"} = $buff;
				$tmp{"name"} = $name;
				$tmp{"desc"} = $desc;
				push(@ret,\%tmp);
				$buff = "";
			}
			$name = $tname;
			$desc = $tdesc;
		}else{
			if($ss =~ />/){
				die $ss;
			}
			$buff .= $ss;
		}
	}
	close(IN);
	if(length($buff) > 0){
		my %tmp;
		$buff =~ s/[\s]//g;
		$tmp{"seq"} = $buff;
		$tmp{"name"} = $name;
		$tmp{"desc"} = $desc;
		push(@ret,\%tmp);
		$buff = "";
	}
	return \@ret;
}

sub arg_to_hash{
	my @args = @{$_[0]};
	my %hashh;
	for(my $ii = 0;$ii <= $#args;$ii++){
		if($args[$ii] =~ /^-/){
			if($ii < $#args){
				$hashh{$args[$ii]} = $args[$ii+1];
			}
		}else{
			$hashh{$args[$ii]} = 1;
		}
	}
	return \%hashh;
}


my %arghash = %{arg_to_hash(\@ARGV)};

my $fastafile = $arghash{"-in"};
my $clusterout = $arghash{"-out"};
my $outfas_dir = $arghash{"-outdir"};
my $ident_threshold = $arghash{"-ident"};
my $cov_long = $arghash{"-cov_long"};
my $cov_short = $arghash{"-cov_short"};

if($ident_threshold > 1){
	$ident_threshold /= 100.0;
}
if($cov_long > 1){
	$cov_long /= 100.0;
}
if($cov_short > 1){
	$cov_short /= 100.0;
}

#my $fastafile = "samples/small_fas.fas";
#my $clusterout = "clustered.dat";
#my $outfas_dir = "clusout/";
#my $ident_threshold = 30;
#my $cov_long = 0.5;
#my $cov_short = 0.8;



if(!-d $outfas_dir){
	mkdir($outfas_dir);
}



my @allseqs = @{get_seq($fastafile)};
@allseqs = sort{length(${$b}{"seq"}) <=> length(${$a}{"seq"}) || ${$a}{"name"} cmp ${$b}{"name"}  }@allseqs;

my $tmpfasname = "tmpin.$$.fas";
my $tmpdbname = "tmpin.$$.db";

open(OUT,">$tmpdbname");
my %seqlength;
my %seq_hs;
foreach my $aa(@allseqs){
	print OUT ">".${$aa}{"name"}."\n";
	print OUT ${$aa}{"seq"}."\n";
	$seqlength{${$aa}{"name"}} = length(${$aa}{"seq"});
	$seq_hs{${$aa}{"name"}} = $aa;
}
close(OUT);
system("$makeblastdbpath -in $tmpdbname -dbtype prot -parse_seqids ");
my $numseqs = $#allseqs+1;
my %clustered;
my $ci = 0;

open(COUT,"> $clusterout");
open(CFOUT,"> $clusterout.fas");
foreach my $qfas(@allseqs){
	if(${$qfas}{"name"} ne "1bb1_B"){
	
		if(defined $clustered{${$qfas}{"name"}}){
			next;
		}
	}
	$ci += 1;
	open(OUT,">$tmpfasname");
	print OUT ">".${$qfas}{"name"}."\n";
	print OUT ${$qfas}{"seq"}."\n";
	close(OUT);
	my $qname = ${$qfas}{"name"};
	#$clustered{$qname}} = $qname;
	
	my $qfilename = $qname;
	$qfilename =~ s/[^A-Za-z0-9\\-\\.]/_/g;
	my $outfas = $outfas_dir."/".$qfilename.".".$ci.".fas";
	
	
	my @ares = `$blastppath -seg no -comp_based_stats 0 -max_target_seqs $numseqs -evalue 100 -outfmt "6 qseqid sseqid evalue qstart qend sstart send pident" -query $tmpfasname -db $tmpdbname `;
	
	my %covered_res;
	my %covered_q;
	
	open(FOUT,"> $outfas.res");
	foreach my $aa(@ares){
		print FOUT $aa;
	}
	close(FOUT);
	foreach my $ress(@ares){
		$ress =~ s/[\r\n]//g;
		my @ptt = split(/\t/,$ress);
		my $qname = $ptt[0];
		my $sname = $ptt[1];
		my $eval = $ptt[2];
		my $qstart = $ptt[3];
		my $qend = $ptt[4];
		my $sstart = $ptt[5];
		my $send = $ptt[6];
		my $pident = $ptt[7];
		
		if($qname =~ /^pdb\|([0-9][0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z])\|(.+)/){
			my $h1 = $1;
			my $h2 = $2;
			$qname = $h1."_".$h2;
		}
		
		if($sname =~ /^pdb\|([0-9][0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z])\|(.+)/){
			my $h1 = $1;
			my $h2 = $2;
			$sname = $h1."_".$h2;
		}
		
		if($pident/100.0 < $ident_threshold){
			next;
		}
		if(!defined $covered_res{$sname}){
			my %tmpp;
			$covered_res{$sname} = \%tmpp;
			my %tmpp2;
			$covered_q{$sname} = \%tmpp2;
		}
		for(my $ii = $qstart;$ii <= $qend;$ii++){
			${$covered_q{$sname}}{$ii} = 1;
		}
		for(my $ii = $sstart;$ii <= $send;$ii++){
			${$covered_res{$sname}}{$ii} = 1;
		}
	}
	open(FOUT,"> $outfas");
	print FOUT ">".${$seq_hs{$qname}}{"name"}." ".${$seq_hs{$qname}}{"desc"}."\n".${$seq_hs{$qname}}{"seq"}."\n";
	print FOUT "\n";
	print CFOUT ">".${$seq_hs{$qname}}{"name"}." ".${$seq_hs{$qname}}{"desc"}."\n".${$seq_hs{$qname}}{"seq"}."\n";
	print CFOUT "\n";
	print COUT $qname;
	foreach my $kk(keys %covered_res){
		if($kk eq $qname){
			next;
		}
		if(defined $clustered{$kk}){
			next;
		}
		my @s1 = keys %{$covered_q{$kk}};
		my @s2 = keys %{$covered_res{$kk}};
		if(!defined $seqlength{$qname}){
			print STDERR $qname."???";
		}
		if(!defined $seqlength{$kk}){
			print STDERR $kk."???";
		}
		if($seqlength{$qname} == $seqlength{$kk}){
			if(
			(($#s1)/$seqlength{$qname} >= $cov_long && ($#s2)/$seqlength{$kk} >= $cov_short)
			|| (($#s1)/$seqlength{$qname} >= $cov_short && ($#s2)/$seqlength{$kk} >= $cov_long)
			){
			}else{
				next;
			}
		}else{
		
			if(($#s1)/$seqlength{$qname} < $cov_long){
				next;
			}
			if(($#s2)/$seqlength{$kk} < $cov_short){
				next;
			}
		}
		
		print FOUT ">".${$seq_hs{$kk}}{"name"}." ".${$seq_hs{$kk}}{"desc"}."\n".${$seq_hs{$kk}}{"seq"}."\n";
		print FOUT "\n";
		$clustered{$kk} = $qname;
		print COUT " ".$kk;
	}
	print COUT "\n";
}
close(COUT);
close(CFOUT);
unlink($tmpfasname);
unlink($tmpdbname);
unlink($tmpdbname.".psi");
unlink($tmpdbname.".phr");
unlink($tmpdbname.".pin");
unlink($tmpdbname.".pog");
unlink($tmpdbname.".psd");
unlink($tmpdbname.".psq");
