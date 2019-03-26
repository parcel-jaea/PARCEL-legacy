all : 
	cd MPI    ; make ;cd ../
	cd SPARSE ; make ;cd ../ 
clean :
	cd MPI    ; make clean ;cd ../
	cd SPARSE ; make clean ;cd ../ 
	rm *~