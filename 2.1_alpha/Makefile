all: compile run

clean:
	rm *.mod && rm program

compile:
	gfortran main.f90 -o program -L /home/andr/inm_ras/ani2D-3.1/lib -lrcb2D-3.1 -lview2D-3.1 -laft2D-3.1 -lfem2D-3.1 -lmba2D-3.1 -lilu-3.1 -lblas 

run:
	./program

