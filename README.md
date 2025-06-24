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
- bloom add only x coordinate uint64_t bits64[3] part

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
- bloom check by x coordinate uint64_t bits64[3] part, computing y coordinate only if there is a hit

Timings are relevant to my PC.
Yours might differ in a great way according to your CPU specs.
  
[alexander@alexander-home Point_Search_JLP]$ ./generate_bloom
[23:28:29] P_table generated
[23:28:29] Range Start: 59 bits
[23:28:29] Range End  : 60 bits
[23:28:29] Block Width: 2^30
[23:28:29] Search Pub : 035c6acbbaa2d43d3134499d937692516202b9de802b739b92b051a05aa4729890
[23:28:29] Settings written to file
[23:28:29] Creating bloomfilter images
[00:04:48] Writing bloom1 image to bloom1.bf
[00:04:52] Writing bloom2 image to bloom2.bf
[00:06:00] Elapsed time: (0)hours (37)minutes (30)seconds

[alexander@alexander-home Point_Search_JLP]$ ./point_search
[00:09:56] S_table generated
[00:09:56] Range Start: 59 bits
[00:09:56] Range End  : 60 bits
[00:09:56] Block Width: 2^30
[00:09:56] Search Pub : 035c6acbbaa2d43d3134499d937692516202b9de802b739b92b051a05aa4729890
[00:09:56] Loading Bloomfilter images
[00:10:21] Search in progress...
[00:10:47] BloomFilter Hit bloom2.bf (Odd Point) [Higher Range Half]
[00:10:47] Private key: 1107029469458650963
[00:10:47] Elapsed time: (0)hours (0)minutes (26)seconds

./generate_bloom uses multiple threads to fill in the bloomfilter binary.
to split the space evenly, number of cores needs to be some power of two value.
actual number of processing cores but equal to some power of two value(2,4,8,16,32,64,...)
divided by 2
actual cores = 8  8 / 2 = 4 cores

./point_search is totally scalable and has no such restriction.

</pre>
