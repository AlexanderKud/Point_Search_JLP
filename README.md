<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
Requires C/C++ OpenMP Library to be installed. <a href="https://www.openmp.org">https://www.openmp.org</a>
  
[alexander@alexander-home Point_Search_JLP]$ ./generate_bloom
[02:00:44] P_table generated
[02:00:44] Range Start: 54 bits
[02:00:44] Range End  : 55 bits
[02:00:44] Block Width: 2^26
[02:00:44] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[02:00:44] Settings written to file
[02:00:44] Creating bloom2 image with 4 threads
[02:00:44] Creating bloom1 image with 4 threads
[02:03:13] Writing bloom2 image to bloom2.bf
[02:03:14] Writing bloom1 image to bloom1.bf
[02:03:14] Elapsed time: (0)hours (2)minutes (31)seconds


[alexander@alexander-home Point_Search_JLP]$ ./point_search
[02:03:20] S_table generated
[02:03:20] Range Start: 54 bits
[02:03:20] Range End  : 55 bits
[02:03:20] Block Width: 2^26
[02:03:20] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[02:03:20] Loading Bloomfilter bloom1.bf
[02:03:20] Loading Bloomfilter bloom2.bf
[02:03:20] Search in progress...
[02:03:33] BloomFilter Hit bloom2.bf (Odd Point) [Higher Range Half]
[02:03:33] Privatekey: 29831168849479125
[02:03:33] Elapsed time: (0)hours (0)minutes (12)seconds


Point batch addition with batch inversion under the hood.
Special thanks to NoMachine for his code draft!!!

Special thanks to kTimesG(<a href="https://github.com/kTimesG">https://github.com/kTimesG</a>) for mutex use improvement proposal.

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
