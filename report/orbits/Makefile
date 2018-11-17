all: orb1.pdf orb2.pdf orb3.pdf orb4.pdf orb5.pdf orb6.pdf orb7.pdf orb8.pdf o1.pdf o2.pdf o3.pdf o4.pdf o5.pdf o6.pdf o7.pdf o8.pdf b1.pdf b2.pdf b3.pdf b4.pdf b5.pdf b6.pdf b7.pdf b8.pdf

.PHONY: clean 
%.pdf: %.eps
	epstopdf --nocompress $*.eps
%.eps: %.gnp %.dat
	gnuplot $*.gnp
%.dat: %.m
	octave -q $*.m > $*.dat

clean: 
	rm *.pdf *.eps *.dat
