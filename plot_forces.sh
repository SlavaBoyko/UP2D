#!/bin/bash

pkill -x gnuplot
rm -rf plots
mkdir plots
cp hat_data.txt plots/

cd plots

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
	
	set title "velosity"	
	set term wxt 2
	set terminal pdf
	set output "velocity.pdf"
 	plot "hat_data.txt" u 1:7 title 'y' with line, \
        "hat_data.txt" u 1:8 title 'x' lt 3 with line

	set title "acceler"	
	set term wxt 2
	set terminal pdf
	set output "acceler.pdf"
 	plot "hat_data.txt" u 1:9 title 'y' with line, \
        "hat_data.txt" u 1:10 title 'x' lt 3 with line

	set title "ang_pos"	
	set term wxt 2
	set terminal pdf
	set output "ang_pos.pdf"
 	plot "hat_data.txt" u 1:11 title 'y' with line
	
	set title "ang_velo"	
	set term wxt 2
	set terminal pdf
	set output "ang_velo.pdf"
 	plot "hat_data.txt" u 1:12 title 'y' with line

	set title "ang_acceler"	
	set term wxt 2
	set terminal pdf
	set output "ang_acceler.pdf"
 	plot "hat_data.txt" u 1:13 title 'y' with line


EOFMarker
