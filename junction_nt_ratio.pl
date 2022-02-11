#!/usr/bin/perl
my %seen_first;
my %seen_last;
my $count;
while (<>){
    next if /^>/;
    chomp;
    my $first = substr($_, 0, 2);
    my $last = substr($_, -2);
    $seen_first{$first}++;
    $seen_last{$last}++;
    $i++;
}


foreach my $element (keys %seen_first){
    print "$element\t", $seen_first{$element}/$i,"\n";
}
print "\n\n\n";
foreach my $element (keys %seen_last){
    print "$element\t", $seen_last{$element}/$i,"\n";
}
