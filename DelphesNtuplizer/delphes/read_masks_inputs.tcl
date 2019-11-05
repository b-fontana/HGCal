#!/usr/bin/tclsh

set NCOLS 10
set NROWS 27
set params { 
    {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99}
    {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99}
    {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99}
    {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99}
    {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99}
    {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99}
    {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99}
    {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99}
    {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99} {-99 -99 -99 -99 -99 -99 -99 -99 -99 -99}
}

set fp [open "../../AnalysisCode/final.csv" r]
gets $fp data
set row_counter 0
while { [gets $fp data] >= 0 } {
    if {$row_counter >= $NROWS} {
	puts "This is iterating too much..."
    }
    puts $data
    set col_counter 0
    set index_split_after 0
    set index_split_before 0
    while { $col_counter < $NCOLS } {
	set index_split_after [string first "," $data [expr $index_split_after + 1]]
	if {$index_split_after == -1} { #last iteration
	    set index_split_after end
	    set var [string range $data [expr $index_split_before + 1] $index_split_after]
	    puts $var
	    lset params $row_counter $col_counter $var
	    set index_split_before $index_split_after
	    break
	}

	if {$col_counter == 0} {
	    set var [string range $data $index_split_before  [expr $index_split_after - 1]]
	} else {
	    set var [string range $data [expr $index_split_before + 1] [expr $index_split_after - 1]]
	}

	puts $var;
	lset params $row_counter $col_counter $var
	set index_split_before $index_split_after
	set col_counter [expr $col_counter + 1] 
    }
    set row_counter [expr $row_counter + 1] 
}
close $fp

puts [lindex $params 0]
puts [lindex $params 1 1]
