$ratio=0.5;
for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i] eq "-m"){
		$MATRIX=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
	elsif($ARGV[$i] eq "-p"){
                $ratio=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
        elsif($ARGV[$i] eq "-o"){
                $out=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
}

if(@ARGV==0){
print "simulate.pl -m \"matrix\" -p \"proportion(p) of coding sequence in genome(0.01<p<0.5)\" -o \"output file\"\n";
exit 1;
}
for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i]!=-1){
		print "simulate.pl -m \"matrix\" -p \"proportion(p) of coding sequence in genome(0.01<p<0.5)\" -o \"output file\"\n";
		exit 1;
        }
}
if($out eq ""){
	print "please put name of output\n";
	exit 1;
}
unless(-f $MATRIX){
        print "can not find file (matrix): $MATRIX\n";
        exit 1;
}
if($ratio<=0 or $ratio>0.5){
        print "the proportion(p) of coding sequence is not correct range. Please use the range(0.01<p<0.5)\n";
        exit 1;
}
if(-f $out){
        unlink $out;
}

$ORDER =5;

%aalist = qw(
atg	M	tat     Y	gca     A       gcc     A
gcg     A       gct     A	tgc     C       tgt     C
gac     D       gat     D	gaa     E       gag     E
ttt     F       ttc     F	gga     G       ggc     G
ggg     G       ggt     G	cac     H       cat     H
ata     I       atc     I	att     I       aaa     K
aag     K       tta     L	ttg     L       cta     L
ctc     L       ctg     L	ctt     L       tgg     W
aac     N       aat     N	cca     P       ccc     P
ccg     P       cct     P	caa     Q       cag     Q
aga     R       agg     R	cga     R       cgc     R
cgg     R       cgt     R	agc     S       agt     S
tca     S       tcc     S	tcg     S       tct     S
aca     T       acc     T	acg     T       act     T
gta     V       gtc     V	gtg     V       gtt     V
tac     Y	taa	.	tag     .       tga     .
);

##
open(IN, $MATRIX);
while(<IN>){
	if(/^Order\s\=\s(\d+),\sFrame\s+\=\s+(\d+)/){
		$order=$1;
		$frame=$2;
	}
	elsif(/^([ATGC]+)\s+\d+\s+\((\d+\.\d{4})/){
		$matrix=$order.$frame;
		$$matrix{$1}=$2;
		push(@orders, $order);
		push(@frames, $frame);
		push(@NTides, $1);
		push(@Percen, $2);
	}
}
close IN;

for(my $i=0; $i<@Percen; $i++){
	if($orders[$i]==$ORDER and $frames[$i]==0){
		$seed_matrix{$NTides[$i]}=$Percen[$i]*10000;
	}
	if($orders[$i]==$ORDER and $frames[$i]==1){
		$matrix_one{$NTides[$i]}=$Percen[$i]*10000;
	}
	if($orders[$i]==$ORDER and $frames[$i]==2){
		$matrix_two{$NTides[$i]}=$Percen[$i]*10000;
	}
	if($orders[$i]==$ORDER and $frames[$i]==3){
		$matrix_thr{$NTides[$i]}=$Percen[$i]*10000;
	}
}

open(out, "> $out");
for(my $len=30; $len<304; $len=$len+3){
	$stop_signal=0;
	undef @ncds_p;
	undef @cds_p;
	for(my $j=0; $j<1000; $j++){
		$num=$j+1;
		$ncds_seq=&random_NCDS_seq($len);
		$t=&PP_calc($ncds_seq);
		push(@ncds_p,$t);

		$cds_seq=&random_CDS_seq($len);
		$t=&PP_calc($cds_seq);
		push(@cds_p,$t);
	}
	@s_temp_p=sort {$b<=>$a} @ncds_p;
	print out "$len\tNCDS\t@s_temp_p\n";
	$SUM=0;
	for(my $n=5; $n<15; $n++){
		$SUM=$SUM+$s_temp_p[$n];
	}
	$AV=$SUM/10;
	$FasleP{$len}=$AV;

	@s_cds_p=sort {$a<=>$b} @cds_p;
	print out "$len\tCDS\t@s_cds_p\n";
	$SUM=0;
	for(my $n=5; $n<15; $n++){
		$SUM=$SUM+$s_temp_p[$n];
	}
	$AV=$SUM/10;
	$FalseN{$len}=$AV;

	$PP=join(",",@s_cds_p);
	$PPcds{$len}=$PP;
}
close out;

sub random_NCDS_seq{
	my ($len, $seq, $order_seq, $next_seq);
	$ref_len=$_[0]; $len=0; undef $seq; undef $order_seq; undef $next_seq;
	while(1){
		$seq=&first_seq($ORDER-1, 0);
		if($seq eq "TAA" or $seq eq "TAG" or $seq eq "TGA"){
		}
		else{
			last;
		}
	}
	while($len < $ref_len){
		$order_seq=substr($seq, $ORDER*-1);
		$next_seq=&random_NCDS_NT($order_seq);
		$seq=$seq.$next_seq;
		if($len%3==0){
			$codon=substr($seq, -3);
			if($codon eq "TAA" or $codon eq "TAG" or $codon eq "TGA"){
				$seq=substr($seq, 0, $len-1);
			}
		}
		$len=length($seq);
			
	}
	return($seq);
}
	

sub random_NCDS_NT{
	my($seed, $seed_seq, $order, $frame, @NT, $seed_next_seq, @temp_freq,  @sum_freq);
	undef @NT; undef $seed_next_seq; undef @sum_freq; undef @temp_freq;
	$seed_seq=$_[0];
	@NT=("A", "T", "G", "C");
	for(my $j=0; $j<@NT; $j++){
		$seed_next_seq=$seed_seq.$NT[$j];
		push(@temp_freq, $seed_matrix{$seed_next_seq});
	}
	$sum=0;
	push(@sum_freq, $sum);
	for(my $i=0; $i<@temp_freq; $i++){
		$sum=$sum+$temp_freq[$i];
		push(@sum_freq, $sum);
	}
	$seed=int(rand($sum))+1;
	$signal=1;
	for(my $i=0; $i<@sum_freq-1; $i++){
		if($sum_freq[$i+1] >= $seed and $seed > $sum_freq[$i]){
			return($NT[$i]);
		}
	}
}
		
sub first_seq{
	my($order, $frame, @temp_freq, @temp_NTides, @sum_freq);
	$order=$_[0];
	$frame=$_[1];
	undef @temp_freq; undef @temp_NTides; undef @sum_freq; 
	for(my $i=0; $i<@Percen; $i++){
		if($orders[$i]==$order and $frames[$i]==$frame){
			if($NTides[$i] ne "TAA" or $NTides[$i] ne "TGA" or $NTides[$i] ne "TAG"){ 
				push(@temp_freq, $Percen[$i]*10000);
				push(@temp_NTides, $NTides[$i]);
			}
		}
	}
	$sum=0;
	push(@sum_freq, $sum);
	for(my $i=0; $i<@temp_freq; $i++){
		$sum=$sum+$temp_freq[$i];
		push(@sum_freq, $sum);
	}
	$seed=int(rand($sum))+1;
	$signal=1;
	for(my $i=0; $i<@sum_freq-1; $i++){
		if($sum_freq[$i+1] >= $seed and $seed > $sum_freq[$i]){
			return($temp_NTides[$i]);
		}
	}
}
		

sub PP_calc{
	my($TARGET, $sign, $first_F, $order, @PFC, $ci);
	$sign=0; undef @PFC;

	$TARGET=shift;

	$sign=&chech_seq($TARGET);
	if($sign==1){
		$ci=0;
		$S=0;
		$target_win=substr($TARGET,$S,30);
		while(length($target_win)>=30){
			for(my $frame=0; $frame<=6; $frame++){
				$first_F=substr($target_win, 0, $ORDER);
				$order=$ORDER-1;
				$matrix=$order.$frame;
				$PFC[$frame]=$$matrix{$first_F};
			}
			for(my $frame=0; $frame<=6; $frame++){
				for(my $k=0; $k<30-$ORDER; $k++){
					$F=substr($target_win, $k, $ORDER);
					$next_NT=substr($target_win, $k+$ORDER, 1);

					if($frame==0){$nextframe=0;}
					elsif(1<=$frame and $frame <=3){$nextframe=($frame-1+$k)%3+1;}
					elsif(4<=$frame and $frame <=6){$nextframe=($frame-4+$k)%3+4;}
					$matrix=$ORDER.$nextframe; $SEG=$F."A";  $probA=$$matrix{$SEG};
					$matrix=$ORDER.$nextframe; $SEG=$F."T";  $probT=$$matrix{$SEG};
					$matrix=$ORDER.$nextframe; $SEG=$F."G";  $probG=$$matrix{$SEG};
					$matrix=$ORDER.$nextframe; $SEG=$F."C";  $probC=$$matrix{$SEG};
	
					if($probA>=0 and $probT>=0 and $probG>=0 and $probC>=0){
						$sigma=$probA+$probT+$probG+$probC;
						if($sigma==0){
							$condprob=0.001;
						}
						elsif($next_NT eq "A"){
							$condprob=$probA/$sigma;
						}
						elsif($next_NT eq "T"){
							$condprob=$probT/$sigma;
						}
						elsif($next_NT eq "G"){
							$condprob=$probG/$sigma;
						}
						elsif($next_NT eq "C"){
							$condprob=$probC/$sigma;
						}
						else{
							die "stop 1\n";
						}
					}		
					else{
						die "stop 2\n";
					}
					$PFC[$frame]=$PFC[$frame]*$condprob;
				}
			}
			for(my $frame=0; $frame<=6; $frame++){
				$A=1/(1-$ratio);
				$B=6/$ratio;
				if($frame==0){ $div = $PFC[$frame]/$A;}
				else{ $div = $div+$PFC[$frame]/$B;}
			}
			@PP=($PFC[0]/$A/$div,$PFC[1]/$B/$div,$PFC[2]/$B/$div,$PFC[3]/$B/$div,$PFC[4]/$B/$div,$PFC[5]/$B/$div,$PFC[6]/$B/$div);
			$ci=$ci+$PFC[1]/$B/$div;
			$S=$S+3;
			$target_win=substr($TARGET,$S,30);
		}
		return($ci);
	}
	else{
		return(-1);
	}
}

sub chech_seq{
        my ($seq, @nu);
        $seq=shift;
        @nu=split(//, $seq);
        $sign=1;
        for(my $i=0; $i<@nu; $i++){
                if($nu[$i]=~/[ATGC]/){
                }
                else{
                        $sign=-1;
                }
        }
        return($sign);
}
			
sub random_CDS_seq{
	my ($len, $seq, $order_seq, $next_seq);
	$ref_len=$_[0]; $len=0; undef $seq; undef $order_seq; undef $next_seq;
	while(1){
		$seq=&first_seq($ORDER-1, 1);
		$codon=substr($seq, 0, 3);
		if($codon eq "TAA" or $codon eq "TAG" or $codon eq "TGA"){
			$seq=substr($seq, 0, $len-1);
		}
		else{
			last;
		}
	}
	while($len < $ref_len){
		$order_seq=substr($seq, $ORDER*-1);
		$next_seq=&random_CDS_NT($order_seq, (length($seq)-($ORDER-3))%3+1);
		$seq=$seq.$next_seq;
		if($len%3==0){
			$codon=substr($seq, -3);
			if($codon eq "TAA" or $codon eq "TAG" or $codon eq "TGA"){
				$seq=substr($seq, 0, $len-1);
			}
		}
		$len=length($seq);
			
	}
	return($seq);
}
	

sub random_CDS_NT{
	my($seed, $seed_seq, $order, $frame, @NT, $seed_next_seq, @temp_freq,  @sum_freq);
	undef @NT; undef $seed_next_seq; undef @sum_freq; undef @temp_freq;
	$seed_seq=$_[0]; $frame=$_[1];
	@NT=("A", "T", "G", "C");
	
	for(my $j=0; $j<@NT; $j++){
		$seed_next_seq=$seed_seq.$NT[$j];
		if($frame==1){ push(@temp_freq, $matrix_one{$seed_next_seq}); }
		elsif($frame==2){ push(@temp_freq, $matrix_two{$seed_next_seq}); }
		elsif($frame==3){ push(@temp_freq, $matrix_thr{$seed_next_seq}); }
	}
	$sum=0;
	push(@sum_freq, $sum);
	for(my $i=0; $i<@temp_freq; $i++){
		$sum=$sum+$temp_freq[$i];
		push(@sum_freq, $sum);
	}
	$seed=int(rand($sum))+1;
	$signal=1;
	for(my $i=0; $i<@sum_freq-1; $i++){
		if($sum_freq[$i+1] >= $seed and $seed > $sum_freq[$i]){
			return($NT[$i]);
		}
	}
}
		
