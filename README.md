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
[14:48:57] Elapsed time: (0)hours (3)minutes (34)seconds

[alexander@alexander-home Point_Search_CPP2]$ ./generate_bloom_2
[00:05:01] P_table generated
[00:05:01] Range Start: 51 bits
[00:05:01] Range End  : 52 bits
[00:05:01] Block Width: 2^25
[00:05:01] Search Pub : 033195139de0331d7a5cab602c4471f728f2e3fb97ed82f593d49ed30ec3c0ba85
[00:05:01] Settings written to file
[00:05:01] Creating BloomFile1 with 4 threads
[00:05:01] Creating BloomFile2 with 4 threads
[00:06:37] Writing BloomFile2 to bloom2.bf
[00:06:38] Writing BloomFile1 to bloom1.bf
[00:06:38] Elapsed time: (0)hours (1)minutes (37)seconds


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
[14:49:33] Elapsed time: (0)hours (0)minutes (28)seconds

</pre>
