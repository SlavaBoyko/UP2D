#!/bin/bash

pkill -x gnuplot

gnuplot -persist <<-EOFMarker
	#fenstergröse
	set terminal wxt size 1200,600 
	set title "Forces"
	set term wxt 0
	plot "Forces.txt" u 1:3 title 'F_x' with line, \
	     "Forces.txt" u 1:2 title 'F_y' lt 3 with line
	set term wxt 1
 	plot "momentum.txt" u 1:2 title 'momentum' with line
EOFMarker
