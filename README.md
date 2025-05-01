<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
  
[alexander@alexander-home Point_Search_CPP2]$ ./generate_bloom
[23:27:40] P_table generated
[23:27:40] Range Start: 51 bits
[23:27:40] Range End  : 52 bits
[23:27:40] Block Width: 2^25
[23:27:40] Search Pub : 033195139de0331d7a5cab602c4471f728f2e3fb97ed82f593d49ed30ec3c0ba85
[23:27:40] Settings written to file
[23:27:40] Creating BloomFile1
[23:27:40] Creating BloomFile2
[23:35:47] Writing BloomFile2 to bloom2.bf
[23:35:47] Writing BloomFile1 to bloom1.bf
[23:36:08] Elapsed time: (0)hours (8)minutes (28)seconds

[alexander@alexander-home Point_Search_CPP2]$ ./point_search
[23:36:15] S_table generated
[23:36:15] Range Start: 51 bits
[23:36:15] Range End  : 52 bits
[23:36:15] Block Width: 2^25
[23:36:15] Search Pub : 033195139de0331d7a5cab602c4471f728f2e3fb97ed82f593d49ed30ec3c0ba85
[23:36:15] Loading Bloomfilter bloom1.bf
[23:36:35] Loading Bloomfilter bloom2.bf
[23:36:54] Search in progress...
[23:37:32] BloomFilter Hit bloom2.bf (Odd Point) [Lower Range Half]
[23:37:32] Privatekey: 3052923631068611
[23:37:32] Elapsed time: (0)hours (0)minutes (37)seconds
</pre>
