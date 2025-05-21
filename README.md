<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
  
[alexander@alexander-home Point_Search_GMP]$ ./generate_bloom
[16:15:29] P_table generated
[16:15:29] Range Start: 59 bits
[16:15:29] Range End  : 60 bits
[16:15:29] Block Width: 2^30
[16:15:29] Search Pub : 038943e72e428a154da75a6ccd0c32118bdf2bfa54077171a24b7418770d276291
[16:15:29] Settings written to file
[16:15:31] Creating BloomFile1 with 4 threads
[16:15:31] Creating BloomFile2 with 4 threads
[16:48:25] Writing BloomFile1 to bloom1.bf
[16:48:26] Writing BloomFile2 to bloom2.bf
[16:49:38] Elapsed time: (0)hours (34)minutes (9)seconds


[alexander@alexander-home Point_Search_GMP]$ ./point_search
[16:53:31] S_table generated
[16:53:31] Range Start: 59 bits
[16:53:31] Range End  : 60 bits
[16:53:31] Block Width: 2^30
[16:53:31] Search Pub : 038943e72e428a154da75a6ccd0c32118bdf2bfa54077171a24b7418770d276291
[16:53:31] Loading Bloomfilter bloom1.bf
[16:53:34] Loading Bloomfilter bloom2.bf
[16:53:37] Search in progress...
[16:53:49] BloomFilter Hit bloom2.bf (Odd Point) [Lower Range Half]
[16:53:49] Privatekey: 0000000000000000000000000000000000000000000000000bf2c3a1a8fea11f
[16:53:49] Elapsed time: (0)hours (0)minutes (11)seconds


./generate_bloom uses multiple threads to fill in the bloomfilter binary.
now has point batch addition implemented.
SPECIAL THANKS TO NoMachine for his code draft!!!
to split the space evenly number of cores needs to be some power of two value.
actual number of processing cores but equal to some power of two value(2,4,8,16,32,64,...)
divided by 2
actual cores = 8  8 / 2 = 4 cores

./point_search is totally scalable and has no such restriction.
can set any desirable number of cores to use but divided by 2.
because we have two search paths : addition and subtraction.
setting cores beyond hardware concurrency will not yield any additional performance.
TODO: point batch addition with bulk inversion to make things faster.

</pre>
