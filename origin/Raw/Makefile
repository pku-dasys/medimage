flag=-Wall -g -m64 -lm

main: AStar.o A.o c_imp.o minimize_EDGE.o minimize_IMAGE.o main.c header.o
	gcc main.c header.o A.o AStar.o c_imp.o minimize_EDGE.o minimize_IMAGE.o \
		$(flag) -o main

A.o: A.c A.h header.o
	gcc $(flag) -c A.c 

AStar.o: AStar.c AStar.h header.o
	gcc $(flag) -c AStar.c

c_imp.o: c_imp.c c_imp.h header.o
	gcc $(flag) -c c_imp.c

minimize_EDGE.o: minimize_EDGE.c minimize_EDGE.h header.o
	gcc $(flag) -c minimize_EDGE.c

minimize_IMAGE.o: minimize_IMAGE.c minimize_IMAGE.h header.o
	gcc $(flag) -c minimize_IMAGE.c
	
header.o: header.c header.h
	gcc $(flag) -c header.c

clean:
	rm -f \
		main \
		*.o \
