gcc -pedantic -c lift.c -o lift.o
gcc -shared -o lift.so lift.o
