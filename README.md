<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
  
[alexander@alexander-home Point_Search_CPP2]$ ./generate_bloom
[14:43:23] P_table generated
[14:43:23] Range Start: 51 bits
[14:43:23] Range End  : 52 bits
[14:43:23] Block Width: 2^25
[14:43:23] Search Pub : 033195139de0331d7a5cab602c4471f728f2e3fb97ed82f593d49ed30ec3c0ba85
[14:43:23] Settings written to file
[14:43:23] Creating BloomFile2
[14:43:23] Creating BloomFile1
[14:48:54] Writing BloomFile2 to bloom2.bf
[14:48:57] Writing BloomFile1 to bloom1.bf
[14:48:57] Elapsed time: (0)hours (5)minutes (34)seconds

[alexander@alexander-home Point_Search_CPP2]$ ./point_search
[14:49:04] S_table generated
[14:49:04] Range Start: 51 bits
[14:49:04] Range End  : 52 bits
[14:49:04] Block Width: 2^25
[14:49:04] Search Pub : 033195139de0331d7a5cab602c4471f728f2e3fb97ed82f593d49ed30ec3c0ba85
[14:49:04] Loading Bloomfilter bloom1.bf
[14:49:04] Loading Bloomfilter bloom2.bf
[14:49:04] Search in progress...
[14:49:08] BloomFilter Hit bloom2.bf (Odd Point) [Lower Range Half]
[14:49:08] False Positive
[14:49:33] BloomFilter Hit bloom2.bf (Odd Point) [Lower Range Half]
[14:49:33] Privatekey: 3052923631068611
[14:49:33] Elapsed time: (0)hours (0)minutes (29)seconds

</pre>
