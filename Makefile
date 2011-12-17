
wor0ket:
		gcc -o wor0ket -Os wor0ket.c -lglut -lGL

.PHONY : clean
 	clean :
		rm wor0ket wor0ket.o
