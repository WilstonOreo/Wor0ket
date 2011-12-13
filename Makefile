
wor0ket:
		gcc -o wor0ket wor0ket.c -O3 -lm -lglut -lGL

.PHONY : clean
 	clean :
		rm wor0ket wor0ket.o
