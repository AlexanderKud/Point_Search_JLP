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
- bloom add only x coordinate

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
[20:37:21] P_table generated
[20:37:21] Range Start: 57 bits
[20:37:21] Range End  : 58 bits
[20:37:21] Block Width: 2^28
[20:37:21] Search Pub : 02d6fba48770c62dbec6e1f88b100dd4d8d213de06cd451c16a12dacdc52d2703d
[20:37:21] Settings written to file
[20:37:21] Creating bloomfilter images
[20:45:30] Writing bloom1 image to bloom1.bf
[20:45:32] Writing bloom2 image to bloom2.bf
[20:45:33] Elapsed time: (0)hours (8)minutes (12)seconds

[alexander@alexander-home Point_Search_JLP]$ ./point_search
[21:25:01] S_table generated
[21:25:01] Range Start: 57 bits
[21:25:01] Range End  : 58 bits
[21:25:01] Block Width: 2^28
[21:25:01] Search Pub : 02d6fba48770c62dbec6e1f88b100dd4d8d213de06cd451c16a12dacdc52d2703d
[21:25:01] Loading Bloomfilter images
[21:25:03] Search in progress...
[21:25:20] BloomFilter Hit bloom2.bf (Odd Point) [Higher Range Half]
[21:25:20] Private key: 247531681757634005
[21:25:20] Elapsed time: (0)hours (0)minutes (16)seconds


./generate_bloom uses multiple threads to fill in the bloomfilter binary.
to split the space evenly, number of cores needs to be some power of two value.
actual number of processing cores but equal to some power of two value(2,4,8,16,32,64,...)
divided by 2
actual cores = 8  8 / 2 = 4 cores

./point_search is totally scalable and has no such restriction.

</pre>
