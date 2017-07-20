CC=gfortran

all: compile run

clean:
	rm src/*.mod && rm src/*.o && rm program

compile:
	cd src; $(CC) -c diffusion.f90 -o diffusion.o
	cd src; $(CC) -c main.f90 -o main.o 
	cd src; $(CC) main.o diffusion.o -o ../program -L /home/andr/inm_ras/ani2D-3.1/lib -lrcb2D-3.1 -lview2D-3.1 -laft2D-3.1 -lfem2D-3.1 -lmba2D-3.1 -lilu-3.1 -lblas 
run:
	./program

