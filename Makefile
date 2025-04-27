default:
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c secp256k1/Int.cpp -o Int.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c secp256k1/Point.cpp -o Point.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c secp256k1/SECP256K1.cpp -o SECP256K1.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c secp256k1/IntMod.cpp -o IntMod.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c secp256k1/Random.cpp -o Random.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c secp256k1/IntGroup.cpp -o IntGroup.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c base58/Base58.cpp -o Base58.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c bech32/Bech32.cpp -o Bech32.o
	g++ -m64 -mssse3 -Wno-write-strings -O1 -o hash/ripemd160.o -c hash/ripemd160.cpp
	g++ -m64 -mssse3 -Wno-write-strings -O1 -o hash/sha256.o -c hash/sha256.cpp
	g++ -m64 -mssse3 -Wno-write-strings -O1 -o hash/ripemd160_sse.o -c hash/ripemd160_sse.cpp
	g++ -m64 -mssse3 -Wno-write-strings -O1 -o hash/sha256_sse.o -c hash/sha256_sse.cpp
	g++ -m64 -mssse3 -Wno-write-strings -O1 -c bloom/bloom.c -o bloom.o
	#g++ -O1 -c util/util.c -o util.o
	g++ -O1 -Wno-write-strings -c generate_bloom.cpp
	g++ -O1 -Wno-write-strings -c point_search.cpp
	g++ -o generate_bloom generate_bloom.o bloom.o Base58.o Bech32.o SECP256K1.o Int.o IntGroup.o IntMod.o Point.o Random.o hash/ripemd160.o hash/ripemd160_sse.o hash/sha256.o hash/sha256_sse.o -lgmp
	g++ -o point_search point_search.o bloom.o Base58.o Bech32.o SECP256K1.o Int.o IntGroup.o IntMod.o Point.o Random.o hash/ripemd160.o hash/ripemd160_sse.o hash/sha256.o hash/sha256_sse.o -lgmp
