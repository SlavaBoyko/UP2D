#!/bin/bash

pkill -x gnuplot

gnuplot -persist <<-EOFMarker
	#fenstergröse
	 #set datafile commentschars "#"
         #show datafile commentschars
         #unset commentschars

	set terminal wxt size 1200,600 
	set title "Forces"

	set term wxt 0
	set terminal pdf
	set output "Forces.pdf"
	plot "hat_data.txt" u 1:3 title 'F_x' with line, \
	     "hat_data.txt" u 1:2 title 'F_y' lt 3 with line
	
	set title "Momentum"	
	set term wxt 1
	set terminal pdf
	set output "moment.pdf"
 	plot "hat_data.txt" u 1:4 title 'momentum' with line

	set title "position"	
	set term wxt 2
	set terminal pdf
	set output "Position.pdf"
 	plot "hat_data.txt" u 1:6 title 'y' with line, \
             "hat_data.txt" u 1:5 title 'x' lt 3 with line

EOFMarker
