ct3d:ct3d.o main.o tracing.o utility.o
	nvcc $^ -o $@ -arch=sm_52 -lboost_system -lboost_filesystem

ct3d.o:ct3d.cpp ct3d.h
	g++ -O3 -c $< -o $@ -std=c++11

tracing.o:tracing.cpp tracing.h
	g++ -O3 -c $< -o $@ -std=c++11

utility.o:utility.cpp utility.h
	g++ -O3 $< -c -o $@ -std=c++11

main.o:main.cu
	nvcc $< -c -o $@ -arch=sm_52
