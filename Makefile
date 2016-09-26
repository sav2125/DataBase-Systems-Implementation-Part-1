all : main

main : main.c
	gcc main.c -lm -lrt -w -o main

clean :
	rm main
