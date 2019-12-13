all: libvector.a

libvector.a: vector.o
	ar rcs libvector.a vector.o

vector.o: vector.cpp vector.h
	g++ -c vector.cpp

clean:
	rm *.a *.o
