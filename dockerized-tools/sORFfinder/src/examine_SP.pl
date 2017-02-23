for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i] eq "-i"){
                $sORF=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
        elsif($ARGV[$i] eq "-l"){
                $LIST=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
        elsif($ARGV[$i] eq "-o"){
                $OUT=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
}

if(@ARGV==0){
print "examine_SP.pl -i \"sORF sequence\" -l \a file including the files of homology\" -o \"output file\"\n";
exit 1;
}
for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i]!=-1){
print "examine_SP.pl -i \"sORF sequence\" -l \a file including the files of homology\" -o \"output file\"\n";
exit 1;
        }
}
if($OUT eq ""){
	print "please put name of output\n";
	exit 1;
}
unless(-f $sORF){
        print "can not find file: $sORF\n";
        exit 1;
}
unless(-f $LIST){
        print "can not find file: $LIST\n";
        exit 1;
}
open(in, $LIST);
while(<in>){
	if(/^(\S+)/){
		$homo_file=$1;
		unless(-f $homo_file){
        		print "can not find file: $homo_file\n";
        		exit 1;
		}
	}
}
close in;
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

for(my $i=0; $i<@NAME; $i++){
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
        $SEQ_aa{$NAME[$i]}=$aa_seq;
}
close out;

open(in, $LIST);
while(<in>){
	if(/^(\S+)/){
		$homo_file=$1;
		open(in_homo, $homo_file);
		while(<in_homo>){
			if(/^>(\S+)/){
				push(@NAME_Q,$1);
			}
			elsif(/(\w+)/){
				push(@SEQU_D,$1);
			}
		}
		close in_homo;
	}
}
close in;

open(OUT, "> $OUT");
for(my $i=0; $i<@NAME; $i++){
	undef @temp_seq;
	for(my $j=0; $j<@NAME_Q; $j++){
		if($NAME[$i] eq $NAME_Q[$j]){
			push(@temp_seq,$SEQU_D[$j]);
		}
	}
	undef %seen;
	@temp_seq=grep(!$seen{$_}++,@temp_seq);
	undef @orf_seq; undef @data_seq;
	if(@temp_seq==0){
		print OUT ">$NAME[$i]\n\n";
		print OUT "No_homolog\n";
		print OUT "\n";
		next;
	}
	print OUT ">$NAME[$i]\n";
	for(my $j=0; $j<@temp_seq; $j++){
		undef $aa_seq;
		for(my $n=0; $n<length($temp_seq[$j])/3; $n++){
			$cod=substr($temp_seq[$j],$n*3,3);
			$cod=~ tr/A-Z/a-z/;
			if( $cod =~ /[atgc]{3}/){
				$aa_seq = $aa_seq.$aalist{$cod};
			}
			elsif( $cod =~ /[a-z]{3}/){$aa_seq=$aa_seq."X";}
			elsif( $cod eq "---"){$aa_seq=$aa_seq."X";}
		}
		open(out, "> tmp1");
		print out ">sORF\n";
		print out "$SEQ_aa{$NAME[$i]}\n";
		print out ">DATA\n";
		print out "$aa_seq\n";
		close out;
		################################
		`clustalw2 tmp1`; #please make a path to clustalw
		###############################
		undef %ALN;
		open(aaALN,"tmp1.aln");
		while(<aaALN>){
			print OUT;
		        if(/^(\S+)(\s\s\s+)((\w|\-)+)/i){
               			$aa = $3;
                		$ALN{$1} = $ALN{$1}.$aa;
        		}
		}
		close aaALN;

		@aa_sORF =split(//,$ALN{"sORF"});
		undef $t_seq; $num_of_gap=0;
		for(my $m = 0; $m < @aa_sORF; $m++){
			if($aa_sORF[$m] =~ /[A-Z]/i){
				$cut_seq = substr($SEQ{$NAME[$i]},($m-$num_of_gap)*3,3);
				$t_seq = $t_seq.$cut_seq;
			}
			elsif($aa_sORF[$m] eq "-"){	
				$t_seq = $t_seq."---";
				$num_of_gap++;
			}
		}
		$t_seq=~ tr/a-z/A-Z/;
		push(@orf_seq,$t_seq);
		@aa_DATA =split(//,$ALN{"DATA"});
		undef $t_seq; $num_of_gap=0;
		for(my $m = 0; $m < @aa_DATA; $m++){
			if($aa_DATA[$m] =~ /[A-Z]/i){
				$cut_seq = substr($temp_seq[$j],($m-$num_of_gap)*3,3);
				$t_seq = $t_seq.$cut_seq;
			}
			elsif($aa_DATA[$m] eq "-"){	
				$t_seq = $t_seq."---";
				$num_of_gap++;
			}
		}
		$t_seq=~ tr/a-z/A-Z/;
		push(@data_seq,$t_seq);
	}
	$R=&estimate_R;
	undef %syss; undef %nons;
	&estimate_sys_non($R);
	($sys_sub,$non_sub,$sys_site,$non_site)=&estimate_sub;
	if($non_site==0){
		$ka=9999;
	}
	else{
		$ka=$non_sub/$non_site;
	}
	if($sys_site==0){
		$ks=9999;
	}
	else{
		$ks=$sys_sub/$sys_site;
	}
	$A=$sys_sub; $B=$sys_site; $C=$non_sub; $D=$non_site;
	$Up=($A+$B+$C+$D)*(($A*$D)-($B*$C))*(($A*$D)-($B*$C));
	$Bo=($A+$B)*($C+$D)*($A+$C)*($B+$D);
	if($Bo!=0){
                $CHI=$Up/$Bo;
                $P=&Statistics::Distributions::chisqrprob(1,$CHI);
        }
	else{
		$P="ne";
	}
	if($ks!=0){
		$PRESSURE=$ka/$ks;
	}
	else{
		$PRESSURE=9999;
	}
	if($PRESSURE<1 and $P<0.05){
		print OUT "ka $ka\tks $ks\t$PRESSURE\t$P\tYES\n";
	}
	else{
		print OUT "ka $ka\tks $ks\t$PRESSURE\t$P\tNO\n";
	}
}
close OUT;

unlink $FASTA_aa;
unlink $FASTA_blast;
unlink "tmp1";
unlink "tmp1.aln";
unlink "tmp1.dnd";

sub estimate_sys_non{
	my $R=shift;
	my @nu=("a", "t", "g", "c");
	my @ref_codons;
	my ($b_trans,$b_tranv,$u_trans,$u_tranv);
	foreach $ref_codon (keys %aalist){
		$sys=0;
		if($ref_codon ne "tag" or $codon ne "taa" or $codon ne "tga"){
			push(@ref_codons, $ref_codon);
			@bara_codon=split(//, $ref_codon);
			$b_trans=0; $b_tranv=0; $u_trans=0; $u_tranv=0;
			for(my $n=0; $n<@nu; $n++){
				if($bara_codon[0] ne $nu[$n]){
					$temp_codon=$nu[$n].$bara_codon[1].$bara_codon[2];
					if($aalist{$temp_codon} ne "."){
						if   ($bara_codon[0] eq "a" and $nu[$n] eq "g") {$b_trans++;}
						elsif($bara_codon[0] eq "g" and $nu[$n] eq "a") {$b_trans++;}
						elsif($bara_codon[0] eq "c" and $nu[$n] eq "t") {$b_trans++;}
						elsif($bara_codon[0] eq "t" and $nu[$n] eq "c") {$b_trans++;}
						else{$b_tranv++;}
					}
					if($aalist{$temp_codon} eq $aalist{$ref_codon}){
						if   ($bara_codon[0] eq "a" and $nu[$n] eq "g") {$u_trans++;}
						elsif($bara_codon[0] eq "g" and $nu[$n] eq "a") {$u_trans++;}
						elsif($bara_codon[0] eq "c" and $nu[$n] eq "t") {$u_trans++;}
						elsif($bara_codon[0] eq "t" and $nu[$n] eq "c") {$u_trans++;}
						else{$u_tranv++;}
					}
				}
			}
			if   ($b_trans==1 and $b_tranv==2){$sys=$sys+($R/(1+$R))*$u_trans+(0.5/(1+$R))*$u_tranv;}
			elsif($b_trans==1 and $b_tranv==1){$sys=$sys+($R/(0.5+$R))*$u_trans+(0.5/(0.5+$R))*$u_tranv;}
			elsif($b_tranv==2)                {$sys=$sys+0.5*$u_tranv;}
			elsif($b_tranv==1)                {$sys=$sys+$u_tranv;}
			elsif($b_trans==1)                {$sys=$sys+$u_trans;}
	
			$b_trans=0; $b_tranv=0; $u_trans=0; $u_tranv=0;
			for(my $n=0; $n<@nu; $n++){
				if($bara_codon[2] ne $nu[$n]){
					$temp_codon=$bara_codon[0].$bara_codon[1].$nu[$n];
					if($aalist{$temp_codon} ne "."){
						if   ($bara_codon[2] eq "a" and $nu[$n] eq "g") {$b_trans++;}
						elsif($bara_codon[2] eq "g" and $nu[$n] eq "a") {$b_trans++;}
						elsif($bara_codon[2] eq "c" and $nu[$n] eq "t") {$b_trans++;}
						elsif($bara_codon[2] eq "t" and $nu[$n] eq "c") {$b_trans++;}
						else{$b_tranv++;}
					}
					if($aalist{$temp_codon} eq $aalist{$ref_codon}){
						if   ($bara_codon[2] eq "a" and $nu[$n] eq "g") {$u_trans++;}
						elsif($bara_codon[2] eq "g" and $nu[$n] eq "a") {$u_trans++;}
						elsif($bara_codon[2] eq "c" and $nu[$n] eq "t") {$u_trans++;}
						elsif($bara_codon[2] eq "t" and $nu[$n] eq "c") {$u_trans++;}
						else{$u_tranv++;}
					}
				}
			}
			if   ($b_trans==1 and $b_tranv==2){$sys=$sys+($R/(1+$R))*$u_trans+(0.5/(1+$R))*$u_tranv;}
			elsif($b_trans==1 and $b_tranv==1){$sys=$sys+($R/(0.5+$R))*$u_trans+(0.5/(0.5+$R))*$u_tranv;}
			elsif($b_tranv==2)                {$sys=$sys+0.5*$u_tranv;}
			elsif($b_tranv==1)                {$sys=$sys+$u_tranv;}
			elsif($b_trans==1)                {$sys=$sys+$u_trans;}
	
			$syss{$ref_codon}=$sys; $nons{$ref_codon}=3-$sys;
		}
	}
}

sub estimate_sub{
	my $S_sub=0; my $N_sub=0;
	my $sum_s_site=0; my $sum_n_site=0;
	for(my $n=0; $n<@orf_seq; $n++){
		my $S_site=0; my $N_site=0;
		for(my $i=0; $i<length($orf_seq[$n])/3; $i++){
			$codon1=substr($orf_seq[$n], $i*3,3);
			$codon2=substr($data_seq[$n],$i*3,3);
			$codon1=~ tr/A-Z/a-z/; $codon2=~ tr/A-Z/a-z/;
			$doing=1;
			if($codon1 =~ /taa|tag|tga/i or $codon2 =~ /taa|tag|tga/i){
				$doing=0;
			}
			if($codon1 =~ /-/i or $codon2 =~ /-/i){
				$doing=0;
			}
			if($doing == 1){
				$t_sys=($syss{$codon1}+$syss{$codon2})/2;
				$t_non=($nons{$codon1}+$nons{$codon2})/2;
				$S_site=$S_site+($syss{$codon1}+$syss{$codon2})/2;
				$N_site=$N_site+($nons{$codon1}+$nons{$codon2})/2;
				@s_n_sub=&sys_non_sub($codon1, $codon2);
				$S_sub=$S_sub+$s_n_sub[0];
				$N_sub=$N_sub+$s_n_sub[1];
			}
		}
		$sum_s_site=$sum_s_site+$S_site;
		$sum_n_site=$sum_n_site+$N_site;
	}
	return($S_sub,$N_sub,$sum_s_site,$sum_n_site)
}

##################ts/tv ratio##################
sub estimate_R{
	my $ts=0; my $tv=0;
	for(my $n=0; $n<@orf_seq; $n++){
		for(my $i=0; $i<length($orf_seq[$n]); $i++){
			$n_t_one=substr($orf_seq[$n], $i, 1);
			$n_t_two=substr($data_seq[$n], $i, 1);
			if($n_t_one ne $n_t_two){
				if($n_t_one eq "A" and $n_t_two eq "G"){$ts++;}
				elsif($n_t_one eq "T" and $n_t_two eq "C"){$ts++;}
				elsif($n_t_one eq "G" and $n_t_two eq "A"){$ts++;}
				elsif($n_t_one eq "C" and $n_t_two eq "T"){$ts++;}
				elsif($n_t_one eq "A" and $n_t_two eq "C"){$tv++;}
				elsif($n_t_one eq "A" and $n_t_two eq "T"){$tv++;}
				elsif($n_t_one eq "T" and $n_t_two eq "A"){$tv++;}
				elsif($n_t_one eq "T" and $n_t_two eq "G"){$tv++;}
				elsif($n_t_one eq "G" and $n_t_two eq "C"){$tv++;}
				elsif($n_t_one eq "G" and $n_t_two eq "T"){$tv++;}
				elsif($n_t_one eq "C" and $n_t_two eq "A"){$tv++;}
				elsif($n_t_one eq "C" and $n_t_two eq "G"){$tv++;}
			}
		}
	}
	if($tv==0){$R=2;}
	else{$R=$ts/$tv;}

	if($R>10 or $R<1){
		$R=2;
	}
	return($R);
}
############################################


#################sys sub, non sub##############################
sub sys_non_sub{
	my ($codon_one, $codon_two, @codon_one_bara, @codon_two_bara);
	my ($sys_sub, $non_sub, $dif);
	
	$codon_one=shift;
	$codon_two=shift;
	@codon_one_bara=split(//, $codon_one);
	@codon_two_bara=split(//, $codon_two);
	$dif=0; $sys_sub=0; $non_sub=0;
	for(my $a=0; $a<3; $a++){
		if($codon_one_bara[$a] ne $codon_two_bara[$a]){
			$dif++;
		}
	}
	if($dif==1){
		if($aalist{$codon_one} eq $aalist{$codon_two}){$sys_sub++;}
		else{$non_sub++;}
	}
	elsif($dif==2){
		$temp_sys=0; $temp_non=0;$div=0;
		if($codon_one_bara[2] eq $codon_two_bara[2]){
			$one_two=$codon_one_bara[0].$codon_two_bara[1].$codon_one_bara[2];
			if($aalist{$one_two} ne "."){
				$div++;
				if($aalist{$codon_one} eq $aalist{$one_two}){$temp_sys++;}
				else{$temp_non++;}
				if($aalist{$one_two} eq $aalist{$codon_two}){$temp_sys++;}
				else{$temp_non++;}
			}

			$one_two=$codon_two_bara[0].$codon_one_bara[1].$codon_one_bara[2];
			if($aalist{$one_two} ne "."){
				$div++;
				if($aalist{$codon_one} eq $aalist{$one_two}){$temp_sys++;}
				else{$temp_non++;}
				if($aalist{$one_two} eq $aalist{$codon_two}){$temp_sys++;}
				else{$temp_non++;}
			}
		}
		elsif($codon_one_bara[1] eq $codon_two_bara[1]){
			$one_two=$codon_one_bara[0].$codon_one_bara[1].$codon_two_bara[2];
			if($aalist{$one_two} ne "."){
				$div++;
				if($aalist{$codon_one} eq $aalist{$one_two}){$temp_sys++;}
				else{$temp_non++;}
				if($aalist{$one_two} eq $aalist{$codon_two}){$temp_sys++;}
				else{$temp_non++;}
			}

			$one_two=$codon_two_bara[0].$codon_one_bara[1].$codon_one_bara[2];
			if($aalist{$one_two} ne "."){
				$div++;
				if($aalist{$codon_one} eq $aalist{$one_two}){$temp_sys++;}
				else{$temp_non++;}
				if($aalist{$one_two} eq $aalist{$codon_two}){$temp_sys++;}
				else{$temp_non++;}
			}
		}
		elsif($codon_one_bara[0] eq $codon_two_bara[0]){
			$one_two=$codon_one_bara[0].$codon_one_bara[1].$codon_two_bara[2];
			if($aalist{$one_two} ne "."){
				$div++;
				if($aalist{$codon_one} eq $aalist{$one_two}){$temp_sys++;}
				else{$temp_non++;}
				if($aalist{$one_two} eq $aalist{$codon_two}){$temp_sys++;}
				else{$temp_non++;}
			}

			$one_two=$codon_one_bara[0].$codon_two_bara[1].$codon_one_bara[2];
			if($aalist{$one_two} ne "."){
				$div++;
				if($aalist{$codon_one} eq $aalist{$one_two}){$temp_sys++;}
				else{$temp_non++;}
				if($aalist{$one_two} eq $aalist{$codon_two}){$temp_sys++;}
				else{$temp_non++;}
			}
		
		}
		$sys_sub=$sys_sub+($temp_sys/$div);
		$non_sub=$non_sub+($temp_non/$div);
	}
	elsif($dif==3){
		$div=0;
		$temp_sys=0; $temp_non=0;
		$one_mid=$codon_one_bara[0].$codon_one_bara[1].$codon_two_bara[2];	
		$mid_two=$codon_one_bara[0].$codon_two_bara[1].$codon_two_bara[2];
		if($aalist{$one_mid} ne "." and $aalist{$mid_two} ne "."){
			$div++;
			if($aalist{$codon_one} eq $aalist{$one_mid}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$one_mid} eq $aalist{$mid_two}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$mid_two} eq $aalist{$codon_two}){$temp_sys++;}
			else{$temp_non++;}
		}

		$one_mid=$codon_one_bara[0].$codon_one_bara[1].$codon_two_bara[2];	
		$mid_two=$codon_two_bara[0].$codon_one_bara[1].$codon_two_bara[2];
		if($aalist{$one_mid} ne "." and $aalist{$mid_two} ne "."){
			$div++;
			if($aalist{$codon_one} eq $aalist{$one_mid}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$one_mid} eq $aalist{$mid_two}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$mid_two} eq $aalist{$codon_two}){$temp_sys++;}
			else{$temp_non++;}
		}

		$one_mid=$codon_one_bara[0].$codon_two_bara[1].$codon_one_bara[2];	
		$mid_two=$codon_one_bara[0].$codon_two_bara[1].$codon_two_bara[2];
		if($aalist{$one_mid} ne "." and $aalist{$mid_two} ne "."){
			$div++;
			if($aalist{$codon_one} eq $aalist{$one_mid}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$one_mid} eq $aalist{$mid_two}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$mid_two} eq $aalist{$codon_two}){$temp_sys++;}
			else{$temp_non++;}
		}

		$one_mid=$codon_one_bara[0].$codon_two_bara[1].$codon_one_bara[2];	
		$mid_two=$codon_two_bara[0].$codon_two_bara[1].$codon_one_bara[2];
		if($aalist{$one_mid} ne "." and $aalist{$mid_two} ne "."){
			$div++;
			if($aalist{$codon_one} eq $aalist{$one_mid}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$one_mid} eq $aalist{$mid_two}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$mid_two} eq $aalist{$codon_two}){$temp_sys++;}
			else{$temp_non++;}
		}

		$one_mid=$codon_two_bara[0].$codon_one_bara[1].$codon_one_bara[2];	
		$mid_two=$codon_two_bara[0].$codon_one_bara[1].$codon_two_bara[2];
		if($aalist{$one_mid} ne "." and $aalist{$mid_two} ne "."){
			$div++;
			if($aalist{$codon_one} eq $aalist{$one_mid}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$one_mid} eq $aalist{$mid_two}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$mid_two} eq $aalist{$codon_two}){$temp_sys++;}
			else{$temp_non++;}
		}

		$one_mid=$codon_two_bara[0].$codon_one_bara[1].$codon_one_bara[2];	
		$mid_two=$codon_two_bara[0].$codon_two_bara[1].$codon_one_bara[2];
		if($aalist{$one_mid} ne "." and $aalist{$mid_two} ne "."){
			$div++;
			if($aalist{$codon_one} eq $aalist{$one_mid}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$one_mid} eq $aalist{$mid_two}){$temp_sys++;}
			else{$temp_non++;}
			if($aalist{$mid_two} eq $aalist{$codon_two}){$temp_sys++;}
			else{$temp_non++;}
		}
		
		$sys_sub=$sys_sub+$temp_sys/$div;
		$non_sub=$non_sub+$temp_non/$div;
	}
	return($sys_sub, $non_sub);
}


package Statistics::Distributions;

use strict;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
use constant PI => 3.1415926536;
use constant SIGNIFICANT => 5; # number of significant digits to be returned

require Exporter;

@ISA = qw(Exporter AutoLoader);
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.
@EXPORT_OK = qw(chisqrdistr tdistr fdistr udistr uprob chisqrprob tprob fprob);
$VERSION = '1.02';

# Preloaded methods go here.
   
sub chisqrdistr { # Percentage points  X^2(x^2,n)
	my ($n, $p) = @_;
	if ($n <= 0 || abs($n) - abs(int($n)) != 0) {
		die "Invalid n: $n\n"; # degree of freedom
	}
	if ($p <= 0 || $p > 1) {
		die "Invalid p: $p\n"; 
	}
	return precision_string(_subchisqr($n, $p));
}

sub udistr { # Percentage points   N(0,1^2)
	my ($p) = (@_);
	if ($p > 1 || $p <= 0) {
		die "Invalid p: $p\n";
	}
	return precision_string(_subu($p));
}

sub tdistr { # Percentage points   t(x,n)
	my ($n, $p) = @_;
	if ($n <= 0 || abs($n) - abs(int($n)) != 0) {
		die "Invalid n: $n\n";
	}
	if ($p <= 0 || $p >= 1) {
		die "Invalid p: $p\n";
	}
	return precision_string(_subt($n, $p));
}

sub fdistr { # Percentage points  F(x,n1,n2)
	my ($n, $m, $p) = @_;
	if (($n<=0) || ((abs($n)-(abs(int($n))))!=0)) {
		die "Invalid n: $n\n"; # first degree of freedom
	}
	if (($m<=0) || ((abs($m)-(abs(int($m))))!=0)) {
		die "Invalid m: $m\n"; # second degree of freedom
	}
	if (($p<=0) || ($p>1)) {
		die "Invalid p: $p\n";
	}
	return precision_string(_subf($n, $m, $p));
}

sub uprob { # Upper probability   N(0,1^2)
	my ($x) = @_;
	return precision_string(_subuprob($x));
}

sub chisqrprob { # Upper probability   X^2(x^2,n)
	my ($n,$x) = @_;
	if (($n <= 0) || ((abs($n) - (abs(int($n)))) != 0)) {
		die "Invalid n: $n\n"; # degree of freedom
	}
	return precision_string(_subchisqrprob($n, $x));
}

sub tprob { # Upper probability   t(x,n)
	my ($n, $x) = @_;
	if (($n <= 0) || ((abs($n) - abs(int($n))) !=0)) {
		die "Invalid n: $n\n"; # degree of freedom
	}
	return precision_string(_subtprob($n, $x));
}

sub fprob { # Upper probability   F(x,n1,n2)
	my ($n, $m, $x) = @_;
	if (($n<=0) || ((abs($n)-(abs(int($n))))!=0)) {
		die "Invalid n: $n\n"; # first degree of freedom
	}
	if (($m<=0) || ((abs($m)-(abs(int($m))))!=0)) {
		die "Invalid m: $m\n"; # second degree of freedom
	} 
	return precision_string(_subfprob($n, $m, $x));
}


sub _subfprob {
	my ($n, $m, $x) = @_;
	my $p;

	if ($x<=0) {
		$p=1;
	} elsif ($m % 2 == 0) {
		my $z = $m / ($m + $n * $x);
		my $a = 1;
		for (my $i = $m - 2; $i >= 2; $i -= 2) {
			$a = 1 + ($n + $i - 2) / $i * $z * $a;
		}
		$p = 1 - ((1 - $z) ** ($n / 2) * $a);
	} elsif ($n % 2 == 0) {
		my $z = $n * $x / ($m + $n * $x);
		my $a = 1;
		for (my $i = $n - 2; $i >= 2; $i -= 2) {
			$a = 1 + ($m + $i - 2) / $i * $z * $a;
		}
		$p = (1 - $z) ** ($m / 2) * $a;
	} else {
		my $y = atan2(sqrt($n * $x / $m), 1);
		my $z = sin($y) ** 2;
		my $a = ($n == 1) ? 0 : 1;
		for (my $i = $n - 2; $i >= 3; $i -= 2) {
			$a = 1 + ($m + $i - 2) / $i * $z * $a;
		} 
		my $b = PI;
		for (my $i = 2; $i <= $m - 1; $i += 2) {
			$b *= ($i - 1) / $i;
		}
		my $p1 = 2 / $b * sin($y) * cos($y) ** $m * $a;

		$z = cos($y) ** 2;
		$a = ($m == 1) ? 0 : 1;
		for (my $i = $m-2; $i >= 3; $i -= 2) {
			$a = 1 + ($i - 1) / $i * $z * $a;
		}
		$p = max(0, $p1 + 1 - 2 * $y / PI
			- 2 / PI * sin($y) * cos($y) * $a);
	}
	return $p;
}


sub _subchisqrprob {
	my ($n,$x) = @_;
	my $p;

	if ($x <= 0) {
		$p = 1;
	} elsif ($n > 100) {
		$p = _subuprob((($x / $n) ** (1/3)
				- (1 - 2/9/$n)) / sqrt(2/9/$n));
	} elsif ($x > 400) {
		$p = 0;
	} else {   
		my ($a, $i, $i1);
		if (($n % 2) != 0) {
			$p = 2 * _subuprob(sqrt($x));
			$a = sqrt(2/PI) * exp(-$x/2) / sqrt($x);
			$i1 = 1;
		} else {
			$p = $a = exp(-$x/2);
			$i1 = 2;
		}

		for ($i = $i1; $i <= ($n-2); $i += 2) {
			$a *= $x / $i;
			$p += $a;
		}
	}
	return $p;
}

sub _subu {
	my ($p) = @_;
	my $y = -log(4 * $p * (1 - $p));
	my $x = sqrt(
		$y * (1.570796288
		  + $y * (.03706987906
		  	+ $y * (-.8364353589E-3
			  + $y *(-.2250947176E-3
			  	+ $y * (.6841218299E-5
				  + $y * (0.5824238515E-5
					+ $y * (-.104527497E-5
					  + $y * (.8360937017E-7
						+ $y * (-.3231081277E-8
						  + $y * (.3657763036E-10
							+ $y *.6936233982E-12)))))))))));
	$x = -$x if ($p>.5);
	return $x;
}

sub _subuprob {
	my ($x) = @_;
	my $p = 0; # if ($absx > 100)
	my $absx = abs($x);

	if ($absx < 1.9) {
		$p = (1 +
			$absx * (.049867347
			  + $absx * (.0211410061
			  	+ $absx * (.0032776263
				  + $absx * (.0000380036
					+ $absx * (.0000488906
					  + $absx * .000005383)))))) ** -16/2;
	} elsif ($absx <= 100) {
		for (my $i = 18; $i >= 1; $i--) {
			$p = $i / ($absx + $p);
		}
		$p = exp(-.5 * $absx * $absx) 
			/ sqrt(2 * PI) / ($absx + $p);
	}

	$p = 1 - $p if ($x<0);
	return $p;
}

   
sub _subt {
	my ($n, $p) = @_;

	if ($p >= 1 || $p <= 0) {
		die "Invalid p: $p\n";
	}

	if ($p == 0.5) {
		return 0;
	} elsif ($p < 0.5) {
		return - _subt($n, 1 - $p);
	}

	my $u = _subu($p);
	my $u2 = $u ** 2;

	my $a = ($u2 + 1) / 4;
	my $b = ((5 * $u2 + 16) * $u2 + 3) / 96;
	my $c = (((3 * $u2 + 19) * $u2 + 17) * $u2 - 15) / 384;
	my $d = ((((79 * $u2 + 776) * $u2 + 1482) * $u2 - 1920) * $u2 - 945) 
				/ 92160;
	my $e = (((((27 * $u2 + 339) * $u2 + 930) * $u2 - 1782) * $u2 - 765) * $u2
			+ 17955) / 368640;

	my $x = $u * (1 + ($a + ($b + ($c + ($d + $e / $n) / $n) / $n) / $n) / $n);

	if ($n <= log10($p) ** 2 + 3) {
		my $round;
		do { 
			my $p1 = _subtprob($n, $x);
			my $n1 = $n + 1;
			my $delta = ($p1 - $p) 
				/ exp(($n1 * log($n1 / ($n + $x * $x)) 
					+ log($n/$n1/2/PI) - 1 
					+ (1/$n1 - 1/$n) / 6) / 2);
			$x += $delta;
			$round = sprintf("%.".abs(int(log10(abs $x)-4))."f",$delta);
		} while (($x) && ($round != 0));
	}
	return $x;
}

sub _subtprob {
	my ($n, $x) = @_;

	my ($a,$b);
	my $w = atan2($x / sqrt($n), 1);
	my $z = cos($w) ** 2;
	my $y = 1;

	for (my $i = $n-2; $i >= 2; $i -= 2) {
		$y = 1 + ($i-1) / $i * $z * $y;
	} 

	if ($n % 2 == 0) {
		$a = sin($w)/2;
		$b = .5;
	} else {
		$a = ($n == 1) ? 0 : sin($w)*cos($w)/PI;
		$b= .5 + $w/PI;
	}
	return max(0, 1 - $b - $a * $y);
}

sub _subf {
	my ($n, $m, $p) = @_;
	my $x;

	if ($p >= 1 || $p <= 0) {
		die "Invalid p: $p\n";
	}

	if ($p == 1) {
		$x = 0;
	} elsif ($m == 1) {
		$x = 1 / (_subt($n, 0.5 - $p / 2) ** 2);
	} elsif ($n == 1) {
		$x = _subt($m, $p/2) ** 2;
	} elsif ($m == 2) {
		my $u = _subchisqr($m, 1 - $p);
		my $a = $m - 2;
		$x = 1 / ($u / $m * (1 +
			(($u - $a) / 2 +
				(((4 * $u - 11 * $a) * $u + $a * (7 * $m - 10)) / 24 +
					(((2 * $u - 10 * $a) * $u + $a * (17 * $m - 26)) * $u
						- $a * $a * (9 * $m - 6)
					)/48/$n
				)/$n
			)/$n));
	} elsif ($n > $m) {
		$x = 1 / _subf2($m, $n, 1 - $p)
	} else {
		$x = _subf2($n, $m, $p)
	}
	return $x;
}

sub _subf2 {
	my ($n, $m, $p) = @_;
	my $u = _subchisqr($n, $p);
	my $n2 = $n - 2;
	my $x = $u / $n * 
		(1 + 
			(($u - $n2) / 2 + 
				(((4 * $u - 11 * $n2) * $u + $n2 * (7 * $n - 10)) / 24 + 
					(((2 * $u - 10 * $n2) * $u + $n2 * (17 * $n - 26)) * $u 
						- $n2 * $n2 * (9 * $n - 6)) / 48 / $m) / $m) / $m);
	my $delta;
	do {
		my $z = exp(
			(($n+$m) * log(($n+$m) / ($n * $x + $m)) 
				+ ($n - 2) * log($x)
				+ log($n * $m / ($n+$m))
				- log(4 * PI)
				- (1/$n  + 1/$m - 1/($n+$m))/6
			)/2);
		$delta = (_subfprob($n, $m, $x) - $p) / $z;
		$x += $delta;
	} while (abs($delta)>3e-4);
	return $x;
}

sub _subchisqr {
	my ($n, $p) = @_;
	my $x;

	if (($p > 1) || ($p <= 0)) {
		die "Invalid p: $p\n";
	} elsif ($p == 1){
		$x = 0;
	} elsif ($n == 1) {
		$x = _subu($p / 2) ** 2;
	} elsif ($n == 2) {
		$x = -2 * log($p);
	} else {
		my $u = _subu($p);
		my $u2 = $u * $u;

		$x = max(0, $n + sqrt(2 * $n) * $u 
			+ 2/3 * ($u2 - 1)
			+ $u * ($u2 - 7) / 9 / sqrt(2 * $n)
			- 2/405 / $n * ($u2 * (3 *$u2 + 7) - 16));

		if ($n <= 100) {
			my ($x0, $p1, $z);
			do {
				$x0 = $x;
				if ($x < 0) {
					$p1 = 1;
				} elsif ($n>100) {
					$p1 = _subuprob((($x / $n)**(1/3) - (1 - 2/9/$n))
						/ sqrt(2/9/$n));
				} elsif ($x>400) {
					$p1 = 0;
				} else {
					my ($i0, $a);
					if (($n % 2) != 0) {
						$p1 = 2 * _subuprob(sqrt($x));
						$a = sqrt(2/PI) * exp(-$x/2) / sqrt($x);
						$i0 = 1;
					} else {
						$p1 = $a = exp(-$x/2);
						$i0 = 2;
					}

					for (my $i = $i0; $i <= $n-2; $i += 2) {
						$a *= $x / $i;
						$p1 += $a;
					}
				}
				$z = exp((($n-1) * log($x/$n) - log(4*PI*$x) 
					+ $n - $x - 1/$n/6) / 2);
				$x += ($p1 - $p) / $z;
				$x = sprintf("%.5f", $x);
			} while (($n < 31) && (abs($x0 - $x) > 1e-4));
		}
	}
	return $x;
}

sub log10 {
	my $n = shift;
	return log($n) / log(10);
}
 
sub max {
	my $max = shift;
	my $next;
	while (@_) {
		$next = shift;
		$max = $next if ($next > $max);
	}	
	return $max;
}

sub min {
	my $min = shift;
	my $next;
	while (@_) {
		$next = shift;
		$min = $next if ($next < $min);
	}	
	return $min;
}

sub precision {
	my ($x) = @_;
	return abs int(log10(abs $x) - SIGNIFICANT);
}

sub precision_string {
	my ($x) = @_;
	if ($x) {
		return sprintf "%." . precision($x) . "f", $x;
	} else {
		return "0";
	}
}


