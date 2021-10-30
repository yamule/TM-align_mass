use strict;
use warnings;


my %abtoletter=("ALA"=>"A",
	"CYS"=>"C",
	"ASP"=>"D",
	"GLU"=>"E",
	"PHE"=>"F",
	"GLY"=>"G",
	"HIS"=>"H",
	"ILE"=>"I",
	"LYS"=>"K",
	"LEU"=>"L",
	"MET"=>"M",
	"ASN"=>"N",
	"PRO"=>"P",
	"GLN"=>"Q",
	"ARG"=>"R",
	"SER"=>"S",
	"THR"=>"T",
	"VAL"=>"V",
	"TRP"=>"W",
	"TYR"=>"Y");

my $targetdir = $ARGV[0];

opendir(DIR,$targetdir);
my @allfiles = grep(/\.pdb$/,readdir(DIR));
closedir(DIR);
my $scount = 0;
my %checked;
foreach my $aa(@allfiles){
	my %cseq = %{getSeq($targetdir."/".$aa)};
	foreach my $ss(keys %cseq){
		my $seqname_base = $aa."_".$ss;
		#my $seqname_base = "seq";
		my $seqname = $seqname_base;
		my $cou = 0;
		while(defined $checked{$seqname}){
			$cou ++;
			$seqname = $seqname_base.".".$cou;
		}
		$checked{$seqname} = 100;
		print ">".$seqname."  file=".$targetdir."/".$aa." chain=$ss  \n";
		print $cseq{$ss}."\n";
	}
}



# 1 -  6        Record name   "ATOM  "
# 7 - 11        Integer       serial       Atom  serial number.
#13 - 16        Atom          name         Atom name.
#17             Character     altLoc       Alternate location indicator.
#18 - 20        Residue name  resName      Residue name.
#22             Character     chainID      Chain identifier.
#23 - 26        Integer       resSeq       Residue sequence number.
#27             AChar         iCode        Code for insertion of residues.
#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#55 - 60        Real(6.2)     occupancy    Occupancy.
#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#77 - 78        LString(2)    element      Element symbol, right-justified.
#79 - 80        LString(2)    charge       Charge  on the atom.
sub getSeq{
	my $filename = $_[0];
	open(IN,$filename) or die "$filename was not found.\n";
	my @lines = <IN>;
	close(IN);
	my %chains;
	my %printed;
	foreach my $ss(@lines){
		if($ss =~ /^TER/){
			last;
		}
		if($ss =~ /^ATOM/ || $ss =~ /^HETATM/){
			$ss =~ s/[\r\n]//g;
			my @pt = split(//,$ss);
			
			my $chain = $pt[21];
			my $res = join("",@pt[17..19]);
			my $rnum = join("",@pt[22..26]);
			
			if(!defined $chains{$chain}){
				$chains{$chain} = "";
			}
			my $code = "X";
			if(!defined $printed{$rnum}){
				if(defined $abtoletter{$res}){
					$code = $abtoletter{$res};
				}
				$chains{$chain} .= $code;
			}
			$printed{$rnum} = 100;
		}
	}
	return \%chains;

}
sub checkChain{
	my $filename = $_[0];
	open(IN,$filename) or die;
	my @lines = <IN>;
	close(IN);
	my %chains;
	foreach my $ss(@lines){
		if($ss =~ /^ATOM/ || $ss =~ /^HETATM/){
			$ss =~ s/[\r\n]//g;
			my @pt = split(//,$ss);
			
			my $chain = $pt[21];
			my $head = join("",@pt[0..29]);
			my $x  = join("",@pt[30..37]);
			my $y  = join("",@pt[38..45]);
			my $z  = join("",@pt[46..53]);
			my $tail = "";
			if($#pt >= 54){
				$tail =  join("",@pt[54..$#pt]);
			}
			my %rec;
			$rec{"head"} = $head;
			$rec{"x"} = $x;
			$rec{"y"} = $y;
			$rec{"z"} = $z;
			$rec{"tail"} = $tail;
			$chains{$chain} = 100;
		}
		#if($ss =~ /^ENDMDL/){
		#	last;
		#}
		#if($ss =~ /^TER[\s]/){
		#	last;
		#}
	}
	return keys %chains;

}
