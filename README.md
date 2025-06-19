<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
Requires C/C++ OpenMP Library to be installed. <a href="https://www.openmp.org">https://www.openmp.org</a>
You can find the right package manager commands on the Internet for your Linux Distro.

On Windows use WSL. Tested to compile and run. (MSYS2 and Cygwin run unstable).

generate_bloom.cpp
- batch addition
- batch inversion
- calculating just x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)

generate_bloom2.cpp
JLP Batch Reference Logic:
- batch addition/subtraction
- batch inversion
- center of the group
- calculating just x coordinate for the batch
- calculating x,y for the last of the batch entry (used as the next startPoint)

point_search.cpp
- batch addition
- batch inversion
- calculating just x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)
- bloom check by x coordinate, computing y coordinate only if there is a hit

Timings are relevant to my PC.
Yours might differ in a great way according to your CPU specs.
  
[alexander@alexander-home Point_Search_JLP]$ ./generate_bloom
[02:00:44] P_table generated
[02:00:44] Range Start: 54 bits
[02:00:44] Range End  : 55 bits
[02:00:44] Block Width: 2^26
[02:00:44] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[02:00:44] Settings written to file
[02:00:44] Creating bloomfilter images
[02:03:13] Writing bloom2 image to bloom2.bf
[02:03:14] Writing bloom1 image to bloom1.bf
[02:03:14] Elapsed time: (0)hours (2)minutes (8)seconds

[alexander@alexander-home Point_Search_JLP]$ ./point_search
[02:03:20] S_table generated
[02:03:20] Range Start: 54 bits
[02:03:20] Range End  : 55 bits
[02:03:20] Block Width: 2^26
[02:03:20] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[02:03:20] Loading Bloomfilter images
[02:03:20] Search in progress...
[02:03:33] BloomFilter Hit bloom2.bf (Odd Point) [Higher Range Half]
[02:03:33] Privatekey: 29831168849479125
[02:03:33] Elapsed time: (0)hours (0)minutes (9)seconds

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
