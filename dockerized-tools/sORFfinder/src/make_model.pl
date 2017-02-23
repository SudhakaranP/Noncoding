$order=5;
for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i] eq "-c"){
                $CDS=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
        elsif($ARGV[$i] eq "-n"){
                $NCDS=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
        elsif($ARGV[$i] eq "-o"){
                $S=0;
                $out=$ARGV[$i+1];
                $ARGV[$i]=-1; $ARGV[$i+1]=-1;
        }
}
if(@ARGV==0){
	print "make_model.pl -c \"coding sequence\" -n \"noncoding sequence\" -o \"output file\"\n";
	exit 1;
}
for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i]!=-1){
		print "@ARGV\n";
		print "make_model.pl -c \"coding sequence\" -n \"noncoding sequence\" -o \"output file\"\n";
		exit 1;
        }
}
if($out eq ""){
	print "please put name of output file\n";
	exit 1;
}
unless(-f $CDS){
        print "can not find file (coding sequence): $CDS\n";
        exit 1;
}
unless(-f $NCDS){
        print "can not find file (noncoding sequence): $NCDS\n";
        exit 1;
}
unless($order==3 or $order==4 or $order==5){
        print "order should be 3 to 5\n";
        exit 1;
}
if(-f $out){
	unlink $out;
}

open(out, ">> $out");

%stop = qw(taa     .	tag     .       tga     .);

open(in, $CDS);
while(<in>){
	if(/^>(\S+)/){
		for(my $i=0; $i <= length($seq)/3; $i++){
			$cod = substr($seq, $i*3,3);
			$cod =~ tr/A-Z/a-z/;
			if( $cod =~ /[atgc]{3}/ and $stop{$cod} ne "."){
				$cod =~ tr/a-z/A-Z/;
				$c_seq = $c_seq.$cod;
			}
		}
		$ant_c_seq=&reverse_seq($c_seq);
		push(@coding, $c_seq);
		push(@A_coding, $ant_c_seq);
		if(@coding==50000){last;}
		undef $seq; undef $c_seq; undef $ant_c_seq;
        }
        elsif(/^(\w+)/){
                $seq = $seq.$1;
	}
}
close in;

for(my $i=0; $i <= length($seq)/3; $i++){
	$cod = substr($seq, $i*3,3);
	$cod =~ tr/A-Z/a-z/;
	if( $cod =~ /[atgc]{3}/ and $stop{$cod} ne "."){
		$c_seq = $c_seq.$cod;
	}
}
$ant_c_seq=&reverse_seq($c_seq);
push(@coding, $c_seq);
push(@A_coding, $ant_c_seq);
if(@coding==0){
	print "there is no conding sequence in $CDS\n";
	exit 1;
}

open(in, $NCDS);
while(<in>){
	if(/^>(\S+)/){
		if($seq ne ""){
			($nc_seq)=&making_orf($seq);
			@t=@$nc_seq;
			if(@t>0){push(@N_coding, @t);}
			if(@N_coding==50000){last;}
		}
		undef $seq; undef $nc_seq;
        }
        elsif(/^(\w+)/){
                $seq = $seq.$1;
	}
}
($nc_seq)=&making_orf($seq);
@t=@$nc_seq;
if(@t>0){push(@N_coding, @t);}

if(@N_coding==0){
	print "there is no noconding sequence in $NCDS\n";
	exit 1;
}

&markov(\@N_coding,($order-1),0);
&markov(\@coding,  ($order-1),1);
&markov(\@coding,  ($order-1),2);
&markov(\@coding,  ($order-1),3);
&markov(\@A_coding,($order-1),4);
&markov(\@A_coding,($order-1),5);
&markov(\@A_coding,($order-1),6);

&markov(\@N_coding,$order,0);
&markov(\@coding,  $order,1);
&markov(\@coding,  $order,2);
&markov(\@coding,  $order,3);
&markov(\@A_coding,$order,4);
&markov(\@A_coding,$order,5);
&markov(\@A_coding,$order,6);

sub making_orf{
	my (@C_SEQ);
	undef @C_SEQ;
	my $seq=shift;
	$frame_seq=substr($seq, 0); $t=&orf60($frame_seq); @t_seq=@$t;push(@C_SEQ, @t_seq);
	$frame_seq=substr($seq, 1); $t=&orf60($frame_seq); @t_seq=@$t;push(@C_SEQ, @t_seq);
	$frame_seq=substr($seq, 2); $t=&orf60($frame_seq); @t_seq=@$t;push(@C_SEQ, @t_seq);

	$rev_seq=&reverse_seq($seq);
	$frame_seq=substr($rev_seq, 0); $t=&orf60($frame_seq); @t_seq=@$t;push(@C_SEQ, @t_seq);
	$frame_seq=substr($rev_seq, 1); $t=&orf60($frame_seq); @t_seq=@$t;push(@C_SEQ, @t_seq);
	$frame_seq=substr($rev_seq, 2); $t=&orf60($frame_seq); @t_seq=@$t;push(@C_SEQ, @t_seq);
	return(\@C_SEQ);
}

sub orf60{
	my($seq, $c_seq, @ORF_seq);
	$seq=shift; undef @ORF_seq; undef $c_seq;
	undef $c_seq; $num=0;
	for(my $i=0; $i <= length($seq)/3; $i++){
		$cod = substr($seq, $i*3,3);
		$cod =~ tr/A-Z/a-z/;
		$not_det=1;
		if($cod =~ /[atgc]{3}/){
			$not_det=0;
		}
		$stop=0;
		if($stop{$cod} eq "."){
			$stop=1;
		}
		$c_seq = $c_seq.$cod;
		if($not_det==1 or $stop==1){
			 if(length($c_seq)>59){
				$num++;
				$SEQ=substr($c_seq, 0, -3) ; $SEQ =~ tr/a-z/A-Z/;
				push(@ORF_seq, $SEQ);
			}
			undef $c_seq;
		}
	}
	return(\@ORF_seq);
}

sub reverse_seq{
	my ($seq, $real_seq, @nor_seq, @rev_seq);
	undef $real_seq; undef @nor_seq; undef @rev_seq;
	$seq=shift;
	@nor_seq=split(//, $seq);
	@rev_seq = reverse @nor_seq;
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
	return ($real_seq);
}

sub markov{
	my (@SEED, %count, $sum, @count_SEED, $order, $ORDER);
	undef @SEED; undef %count; $sum=0; undef @count_SEED;
	my $S=shift; @seq=@$S;
	my $ORDER=shift; my $frame=shift; 
	@NT=('A', 'T', 'G', 'C');
	$order=$ORDER+1;
	$num=4**$order;
	print out "Order = $ORDER, Frame = $frame\n";
	for(my $i=0; $i<$num; $i++){
		undef $seed; $tmp=$i;
		for(my $n=0; $n<$order; $n++){
			$se=$tmp%4; $seed=$seed.$NT[$se];
			$tmp=int($tmp/4);
		}
		push(@SEED, $seed);
	}	
	if($frame==0){
		for(my $i=0; $i<@seq; $i++){
			for(my $j=0; $j<length($seq[$i]); $j++){
				$temp_seq=substr($seq[$i], $j, $order);
				if(length($temp_seq)<$order){last;}
				$count{$temp_seq}++;
			}
		}
		$sum=0;
		for(my $i=0; $i<@SEED; $i++){
			$sum=$sum+$count{$SEED[$i]};	
		}
		for(my $i=0; $i<@SEED; $i++){
			$rate=100*$count{$SEED[$i]}/$sum;
			printf out "%-8s %13d ", $SEED[$i], $count{$SEED[$i]};
			printf out "\(%1.9f\%\)\n", $rate;
		}
	}
	else{
		if   ($frame==4){$frame=1;}
		elsif($frame==5){$frame=2;}
		elsif($frame==6){$frame=3;}

		for(my $i=0; $i<@seq; $i++){
			for(my $j=0; $j<length($seq[$i]); $j=$j+3){
				$temp_seq=substr($seq[$i], ($j+$frame-1), $order);
				if(length($temp_seq)<$order){last;}
				$count{$temp_seq}++;
			}
		}
		$sum=0;
		for(my $i=0; $i<@SEED; $i++){
			$sum=$sum+$count{$SEED[$i]};	
		}
		for(my $i=0; $i<@SEED; $i++){
			$rate=100*$count{$SEED[$i]}/$sum;
			printf out "%-8s %13d ", $SEED[$i], $count{$SEED[$i]};
			printf out "\(%1.9f\%\)\n", $rate;
		}
	}
}

