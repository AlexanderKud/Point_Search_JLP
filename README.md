<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
  
[alexander@alexander-home Point_Search_JLP]$ ./generate_bloom
[21:48:59] P_table generated
[21:48:59] Range Start: 54 bits
[21:48:59] Range End  : 55 bits
[21:48:59] Block Width: 2^26
[21:48:59] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[21:48:59] Settings written to file
[21:48:59] Creating bloom2 image with 4 threads
[21:48:59] Creating bloom1 image with 4 threads
[21:50:53] Writing bloom1 image to bloom1.bf
[21:50:54] Writing bloom2 image to bloom2.bf
[21:50:54] Elapsed time: (0)hours (1)minutes (55)seconds


[alexander@alexander-home Point_Search_JLP]$ ./point_search
[21:51:09] S_table generated
[21:51:09] Range Start: 54 bits
[21:51:09] Range End  : 55 bits
[21:51:09] Block Width: 2^26
[21:51:09] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[21:51:09] Loading Bloomfilter bloom1.bf
[21:51:09] Loading Bloomfilter bloom2.bf
[21:51:09] Search in progress...
[21:51:53] BloomFilter Hit bloom2.bf (Odd Point) [Higher Range Half]
[21:51:53] Privatekey: 29831168849479125
[21:51:53] Elapsed time: (0)hours (0)minutes (44)seconds

[alexander@alexander-home Point_Search_JLP]$ ./point_search_batch
[21:51:58] S_table generated
[21:51:58] Range Start: 54 bits
[21:51:58] Range End  : 55 bits
[21:51:58] Block Width: 2^26
[21:51:58] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[21:51:58] Loading Bloomfilter bloom1.bf
[21:51:58] Loading Bloomfilter bloom2.bf
[21:51:59] Search in progress...
[21:52:25] BloomFilter Hit bloom2.bf (Odd Point) [Higher Range Half]
[21:52:25] Privatekey: 29831168849479125
[21:52:25] Elapsed time: (0)hours (0)minutes (26)seconds

Point batch addition with batch inversion under the hood.
SPECIAL THANKS TO NoMachine for his code draft!!!

./generate_bloom uses multiple threads to fill in the bloomfilter binary.
to split the space evenly, number of cores needs to be some power of two value.
actual number of processing cores but equal to some power of two value(2,4,8,16,32,64,...)
divided by 2
actual cores = 8  8 / 2 = 4 cores

./point_search is totally scalable and has no such restriction.
can set any desirable number of cores to use but divided by 2.
because we have two search paths : addition and subtraction.
setting cores beyond hardware concurrency will not yield any additional performance.

</pre>
