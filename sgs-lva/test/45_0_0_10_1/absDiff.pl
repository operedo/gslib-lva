
open(IN1,"<$ARGV[0]");
open(IN2,"<$ARGV[1]");

$tol=1.0E-1;

$line=<IN1>; # Title
$line=<IN1>; # num vars (1)
$line=<IN1>; # name first var (sim)

$line=<IN2>; # Title
$line=<IN2>; # num vars (1)
$line=<IN2>; # name first var (sim)

$counter=0;
while(<IN1>){

	$line1=$_;
	$line2=<IN2>;
	
	$line1=~/([0-9.E+\-]+)/;
	$est1=$1;
	
	$line2=~/([0-9.E+\-]+)/;
	$est2=$1;

	if(abs($est1-$est2)>$tol){
		$counter=$counter+1;
	}	

}
close(IN1);
close(IN2);
if($counter>0){
	print $counter,"\n";
}
