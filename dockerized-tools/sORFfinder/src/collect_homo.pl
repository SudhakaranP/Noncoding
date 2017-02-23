for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i] eq "-d"){
                $database=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
        elsif($ARGV[$i] eq "-i"){
                $sORF=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
        elsif($ARGV[$i] eq "-o"){
                $OUT=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
}

if(@ARGV==0){
	print "collect_homo.pl -d \"blast database(nucleotide)\" -i \"sORF sequence\" -o \"output file\"\n";
	exit 1;
}
for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i]!=-1){
		print "collect_homo.pl -d \"blast database(nucleotide)\" -i \"sORF sequence\" -o \"output file\"\n";
		exit 1;
        }
}
if($OUT eq ""){
	print "please put name of output\n";
	exit 1;
}
unless(-f $sORF){
        print "can not find file (sORF sequence): $sORF\n";
        exit 1;
}
unless(-f $database){
        print "can not find file (blast database): $database\n";
        exit 1;
}
if(-f $OUT){
        unlink $OUT;
}

%aalist = qw(
atg     M       tat     Y       gca     A       gcc     A
gcg     A       gct     A       tgc     C       tgt     C
gac     D       gat     D       gaa     E       gag     E
ttt     F       ttc     F       gga     G       ggc     G
ggg     G       ggt     G       cac     H       cat     H
ata     I       atc     I       att     I       aaa     K
aag     K       tta     L       ttg     L       cta     L
ctc     L       ctg     L       ctt     L       tgg     W
aac     N       aat     N       cca     P       ccc     P
ccg     P       cct     P       caa     Q       cag     Q
aga     R       agg     R       cga     R       cgc     R
cgg     R       cgt     R       agc     S       agt     S
tca     S       tcc     S       tcg     S       tct     S
aca     T       acc     T       acg     T       act     T
gta     V       gtc     V       gtg     V       gtt     V
tac     Y       taa     .       tag     .       tga     .
);

open(in, $sORF);
while(<in>){
	if(/^#/){next;}
	if(/^>(\S+)/){
		$name=$1;
		if($done{$name}==1){
			print "error! use redundant seq name\n";
			exit 1;
		}
		push(@NAME,$name);
		$done{$name}=1;
	}
	elsif(/(\S+)/){
		$SEQ{$name}=$SEQ{$name}.$1;
	}
}
close in;

$FASTA_aa="tmp.aa";
open(out, "> $FASTA_aa");
for(my $i=0; $i<@NAME; $i++){
	print out ">$NAME[$i]\n";
	undef $aa_seq;
	for(my $j=0; $j<length($SEQ{$NAME[$i]})/3; $j++){
		$cod=substr($SEQ{$NAME[$i]},$j*3,3);
		$cod=~ tr/A-Z/a-z/;
		if( $cod =~ /[atgc]{3}/){
			$aa_seq = $aa_seq.$aalist{$cod};
		}
		elsif( $cod =~ /[a-z]{3}/){
			print "error! some nucleotides does not have A,T,G,C in sORF sequences\n";
		}
		elsif( $cod eq "---"){
			print "error! some nucleotides does not have A,T,G,C in sORF sequences\n";
		}
	}
	print out "$aa_seq\n";
	$SEQ_aa{$NAME[$i]}=$aa_seq;
}
close out;

$FASTA_blast="tmp.blast";
		################################
`tblastn -query $FASTA_aa -db $database -num_alignments 5 -outfmt 6 -out $FASTA_blast`;
		###############################
		#AT1G01230-AT1G01240#98#+#1.3    chr05   35.29   34      22      0       39      72      10815546        10815647        5.9     27.7
		#AT1G01700-AT1G01710#596#+#51.2  chr10   50.00   20      10      0       2       21      4355211 4355152 1.4     28.5
		#AT1G02080-AT1G02090#205#+#53.9  chr01   71.43   14      4       0       2       15      32624729        32624770        0.35    30.4
		#AT1G02190-AT1G02205#595#+#16.6  chr03   32.35   34      23      0       1       34      12681686        12681585        0.54    30.4
		#AT1G02210-AT1G02220#137#+#43.4  chr02   47.83   23      12      0       7       29      13567679        13567747        0.19    31.2
		#AT1G02430-AT1G02440#369#+#82.4  chr03   40.91   44      26      0       1       44      5698664 5698795 6e-06   45.8

open(in, $FASTA_blast);
while(<in>){
	chomp;
	@tmp=split(/\s+/,$_);
	push(@nameQ, $tmp[0]); push(@fiveQ,  $tmp[6]); push(@threQ,  $tmp[7]);
	push(@nameD, $tmp[1]); push(@fiveD,  $tmp[8]); push(@threD,  $tmp[9]);
}
close in;

open(out, "> $OUT");
for(my $i=0; $i<@nameD; $i++){
	if($count{$nameD[$i]}>1){next;}
	$count{$nameD[$i]}++;
	if($fiveD[$i] < $threD[$i]){
		$seq_name=$nameD[$i]."|+|".$fiveD[$i]."-".$threD[$i];
		$GEN_SEQ=&get_genome($nameD[$i], $fiveD[$i], $threD[$i]);
	}
	else{
		$seq_name=$nameD[$i]."|-|".$threD[$i]."-".$fiveD[$i];
		$temp=&get_genome($nameD[$i], $threD[$i], $fiveD[$i]);
		$GEN_SEQ=&rev_seq($temp);
	}
	if($SEQ{$nameQ[$i]}=~/$GEN_SEQ/){next;}
	print out ">$nameQ[$i]\n";
	print out "$GEN_SEQ\n";
}

sub get_genome{
	my ($seq, $len);
	undef $seq; undef $len;
	my $chr=shift;
	my $L=shift;
	my $R=shift;

	$sign=0; $len=0;
	open(in, $database);
	while(<in>){
		if(/>$chr/){
			$sign=1;
			$start=0;
		}
		elsif($sign==1 and /(\w+)/){
			$len=$len+length($1);
			$temp=$1;
			if($len>=$L and $start==0){ $start=1; $startL=$len-length($temp);}
			if($start==1){ $seq=$seq.$temp;}
			if($len>=$R){$endL=$len; last;}
		}
	}
	close in;

	$start=$L-$startL-1;
	$len_s=$R-$L+1;
	$seg=substr($seq, $start, $len_s);
	return($seg);
}

sub rev_seq{
	my($seq, @sequence, @rev_seq, $real_seq);
	$seq=shift; undef $real_seq;
	@sequence=split("",$seq);
	@rev_seq = reverse @sequence;
	for(my $i =0; $i < @rev_seq; $i++){
		if(($rev_seq[$i] eq "A") or ($rev_seq[$i] eq "a")){
			$real_seq = $real_seq."T";	
		}
		if(($rev_seq[$i] eq "T") or ($rev_seq[$i] eq "t")){
			$real_seq = $real_seq."A";	
		}
		if(($rev_seq[$i] eq "G") or ($rev_seq[$i] eq "g")){
			$real_seq = $real_seq."C";	
		}
		if(($rev_seq[$i] eq "C") or ($rev_seq[$i] eq "c")){
			$real_seq = $real_seq."G";	
		}
	}
	return($real_seq);
}

