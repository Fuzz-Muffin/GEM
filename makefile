e = toms707.o gnufor2.o para.o input_class.o grid_class.o basis_class.o basis_class_laguerre.o HTVS_class.o eigen_class.o oscillator_strength.o main.o

a = -C -g -ffpe-trap='underflow'
b = -C -std=legacy

main : $e
	gfortran $a $e -o main -l lapack

gnufor2.o : gnufor2.f90
	gfortran $a -c gnufor2.f90

toms707.o : toms707.f90
	gfortran $a -c toms707.f90

para.o : para.f90 toms707.f90
	gfortran $a -c para.f90

input_class.o : input_class.f90 para.f90
	gfortran $a -c input_class.f90

grid_class.o : grid_class.f90 input_class.f90
	gfortran $a -c grid_class.f90

basis_class.o : basis_class.f90 grid_class.f90 para.f90
	gfortran $a -c basis_class.f90

basis_class_laguerre.o : basis_class_laguerre.f90 grid_class.f90 para.f90 
	gfortran $a -c basis_class_laguerre.f90

HTVS_class.o : HTVS_class.f90 para.f90 basis_class.f90
	gfortran $a -c HTVS_class.f90

eigen_class.o : eigen_class.f90 para.f90
	gfortran $a -c eigen_class.f90

oscillator_strength.o : oscillator_strength.f90 para.f90 grid_class.f90 basis_class.f90 basis_class_laguerre.f90 eigen_class.f90
	gfortran $a -c oscillator_strength.f90

main.o : main.f90 oscillator_strength.f90 eigen_class.f90 HTVS_class.f90 basis_class.f90 basis_class_laguerre.f90 input_class.f90 para.f90 grid_class.f90 gnufor2.f90 toms707.f90
	gfortran $a -c main.f90 -l lapack

clean :
	rm -f *.o *.mod
