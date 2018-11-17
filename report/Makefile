all:potential.pdf 

.PHONY: clean

potential.pdf: potential.tex 
	$(MAKE) -C orbits
	pdflatex potential

clean: 
	rm *.aux *.log *~
