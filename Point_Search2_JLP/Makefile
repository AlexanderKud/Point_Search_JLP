default:
	g++ -m64 -mssse3 -Wno-write-strings -O2 -c secp256k1/Int.cpp -o Int.o
	g++ -m64 -mssse3 -Wno-write-strings -O2 -c secp256k1/Point.cpp -o Point.o
	g++ -m64 -mssse3 -Wno-write-strings -O2 -c secp256k1/SECP256K1.cpp -o SECP256K1.o
	g++ -m64 -mssse3 -Wno-write-strings -O2 -c secp256k1/IntMod.cpp -o IntMod.o
	g++ -m64 -mssse3 -Wno-write-strings -O2 -c secp256k1/IntGroup.cpp -o IntGroup.o
	g++ -O2 -march=native -c util/util.cpp -o util.o
	g++ -std=c++20  -O2 -march=native -mtune=native -funroll-loops -fomit-frame-pointer -Wno-write-strings -c bloom.cpp
	g++ -std=c++20 -O2 -march=native -mtune=native -funroll-loops -fomit-frame-pointer -Wno-write-strings -c search.cpp
	g++ -o bloom bloom.o util.o SECP256K1.o Int.o IntGroup.o IntMod.o Point.o
	g++ -o search search.o util.o SECP256K1.o Int.o IntGroup.o IntMod.o Point.o
	rm *.o
