default:
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c secp256k1/Int.cpp -o Int.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c secp256k1/Point.cpp -o Point.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c secp256k1/SECP256K1.cpp -o SECP256K1.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c secp256k1/IntMod.cpp -o IntMod.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c secp256k1/IntGroup.cpp -o IntGroup.o
	g++ -O1 -Wno-write-strings -c generate_bloom.cpp
	g++ -O1 -Wno-write-strings -c point_search.cpp
	g++ -o generate_bloom generate_bloom.o SECP256K1.o Int.o IntGroup.o IntMod.o Point.o -lgmp
	g++ -o point_search point_search.o SECP256K1.o Int.o IntGroup.o IntMod.o Point.o -lgmp
	rm *.o
