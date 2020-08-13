#!/usr/bin/perl
$infile=$ARGV[0];

(open (INFILE, $infile))||die "ERROR opening $infile";

while($cont=<INFILE>){
    chop($cont);
    if($cont =~ /^Tree length/){
	last;
    }
}
$cont=<INFILE>;

#now start reading the tree;
$i=0;
while($cont=<INFILE>){
    if($cont =~ /Processing/){
	last;
    }
    @{$data[$i]}=split("",$cont);
#    print "line $i:  $data[$i][0] $data[$i][1] $data[$i][2]\n";
    $i++;
}

$num_rows=$i;
#count through rows
for($row=0;$row<$num_rows;$row++){
#    print "now on row $row\n";
    $last=$#{$data[$row]};
    $name="";
    $query="";
    #count through columns from right to left
    for($col=$last-1; $col>=0; $col--){
	$tag1=1;
	$char=$data[$row][$col];
#	print "col is $col, char is $char;\n";
	if($col<1 && $char =~ /[A-Za-z0-9]/){
	    $name.=$char;
#	    print "root touches edge for $name\n"; 
#print output with the correct spacing
	    @b=split("",$name);
	    for($q=$#b;$q>=0;$q--){
		$query.=$b[$q];
	    }
	    $diff=10-length($query);
	    print "$query";
	    for($n=1;$n<=$diff;$n++){
		print " ";
	    }
	    print "\tnone      \n";
#end of print output with correct spacing
	    $name="";
	    $query="";
	    undef $root;
	}
	if($char =~ /[A-Za-z0-9]/){
	    $name.=$char;
	}
	else{
	    if ($name =~ /[A-Za-z0-9]/){
#		print "\n\n***********\n$name is at position $row, $col\n";
		@b=split("",$name);
		for($q=$#b;$q>=0;$q--){
#		    print "$b[$q]";
		    $query.=$b[$q];
		}
#		print "\t";
#		print "$name\n";
#		print "\n\n************************\n**calling track_node for $query at position $col, character '$data[$row][$col]'\n";
		#$col is the position immediately after the name, read from right to left
		$testroot=0;
		for($y=$col;$y>=0;$y--){
		    $testrootchar=$data[$row][$y];
		    if($testrootchar eq "-"){
		    }
		    else{
			$testroot=1;
			last;
		    }
		}
		if($testroot <1){    #is this the root node?
#print output with the correct spacing
		    $diff=10-length($query);
		    print "$query";
		    for($n=1;$n<=$diff;$n++){
			print " ";
		    }
		    print "\tnone      \n";
#end of print output with correct spacing
		}
		else{
		    &track_node;
		}
		$name="";
		$query="";
		undef $root;
	    }
	}
    }
}
exit;

sub track_node{
    #now track this node back to find its ancestor
    $mark_set=0;  #set this=1 once you have found the proper value for $col2 to go up or down
    #first find the point to start looking up or down
    $char=$data[$row][$col];

    if( ($data[$row][$col-1] eq "-") || ($char eq "-")){
#	print "simple non-squished case, look left to find turning point\n";
	#this is the simple non-squished case-look left to find turning point
	for($col2=$col-1; $col2>=0; $col2--){
	    #col2 starts counting at the position immediately after the name, counting right to left 
	    $char=$data[$row][$col2];
#	    print "char is '$char' for $col2\n";
	    if($char =~ /\//){
		$mark_set=1;
		last;
	    }
	    elsif($char =~ /\\/){
		$mark_set=1;
		last;
	    }
	    elsif($char =~ /\+/){
		$mark_set=1;
		last;
	    }
	    elsif($char =~ /[0-9]/ && $data[$row][$col2+1] eq "-"){
		while($char =~ /[0-9]/){
		    $root=$char.$root;
#		    print "growing $root for $query\n";
		    $col2--;
		    $char=$data[$row][$col2];
		}
#print output with the correct spacing
		$diff=10-length($query);
		$diff2=10-length($root);
		if($root){
		    print "$query";
		    for($n=1;$n<=$diff;$n++){
			print " ";
		    }
		    print "\t$root";
		    for($n=1;$n<=$diff2;$n++){
			print " ";
		    }
		    print "\n";
		}
#end of print output with correct spacing

	    }
	}
    }
    elsif($char eq "|" || $char eq " "){
	#this is a very squished case, will need to readjust turning point
#	print "squished case, readjust turning point\n";
	for($col2=$col+1; $col2<$col+5; $col2++){
	    if( ($data[$row+1][$col2]=~/\S/) || ($data[$row-1][$col2] =~ /\S/) ){
		$mark_set=1;
		last;
	    }
	}     
    }    
    elsif($char eq "\\" || $char eq "\/"){
	#this is a regular squished case, will need to readjust turning point
#	print "squished case, readjust turning point\n";
	for($col2=$col; $col2<$col+6; $col2++){
	    if( ($data[$row+1][$col2]=~/\S/) || ($data[$row-1][$col2] =~ /\S/) ){
		$mark_set=1;
		last;
	    }
	}     
    }    
    if($mark_set){
#	print "mark point selected at $col2\n";
#	print "char for $col2 is $data[$row][$col2]\n";
#	print "up:  '$data[$row-1][$col2]' for $row -1, $col2\n";
#	print "down:  '$data[$row+1][$col2]' for $row +1, $col2\n";
	if(($data[$row-1][$col2]=~ /\S/) && ($row > 0)){
	    undef $root;
#	    print "$row, $col2 going up exists:$data[$row-1][$col2]\n";
	    for($row2=$row-1; $row2>=0; $row2--){
		$char=$data[$row2][$col2];
#		print "$char\n";
		if($char =~ /\S/){
		}
		else{
		    last;
		}
		if($char =~ /[0-9a-zA-Z]/){
#		    print "$data[$row2][$col2-2] $data[$row2][$col2-1] $char\n";
#		    print "$char\n";

		    $root_string="";
		    for($a=$col2-5;$a<=$col2+3;$a++){
			$root_string.=$data[$row2][$a];
		    }
		    #works only for internal node numbers in the range 1000-99999
		    if($root_string =~ /(\d{2,5})/){
			$root=$1;
		    }
#print output with the correct spacing
		    $diff=10-length($query);
		    $diff2=10-length($root);
		    if($root){
			print "$query";
			for($n=1;$n<=$diff;$n++){
			    print " ";
			}
			print "\t$root";
			for($n=1;$n<=$diff2;$n++){
			    print " ";
			}
			print "\n";
		    }
#end of print output with correct spacing
		    last;
		}
	    }		
	}
	if(($data[$row+1][$col2]=~ /\S/) && ($row >= 0)){
	    undef $root;
#	    print "going down exists\n";
	    for($row2=$row+1; $row2<=$num_rows; $row2++){
		$char=$data[$row2][$col2];
#		print "$char\n";
		if($char =~ /\S/){
		}
		else{
		    last;
		}
		if($char =~ /[0-9a-zA-Z]/){
#		    print "a$data[$row2][$col2-2] b$data[$row2][$col2-1] c$char d$data[$row2][$col2+1] e$data[$row2][$col2+2]\n";

		    $root_string="";
		    for($a=$col2-5;$a<=$col2+3;$a++){
			$root_string.=$data[$row2][$a];
		    }
		    if($root_string =~ /(\d{2,5})/){
			$root=$1;
		    }
#print the output with correct spacing
		    if($root){
			$diff=10-length($query);
			$diff2=10-length($root);
			print "$query";
			for($n=1;$n<=$diff;$n++){
			    print " ";
			}
			print "\t$root";
			for($n=1;$n<=$diff2;$n++){
			    print " ";
			}
			print "\n";
		    }
#end of print output with correct spacing
		    last;
		}
	    }
	}
    }
}



