<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
Requires C/C++ OpenMP Library to be installed. <a href="https://www.openmp.org">https://www.openmp.org</a>
You can find the right package manager commands on the Internet for your Linux Distro.

On Windows use WSL. Tested to compile and run.

generate_bloom.cpp
- batch addition
- batch inversion
- calculating x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)
- bloom add x coordinate uint64_t bits64[3] part

generate_bloom2.cpp
JLP Batch Reference Logic:
- batch addition/subtraction
- batch inversion
- center of the group
- calculating x coordinate for the batch
- calculating x,y for the last of the batch entry (used as the next startPoint)

point_search.cpp
- batch addition
- batch inversion
- calculating x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)
- bloom check by x coordinate uint64_t bits64[3] part, computing y coordinate only if there is a hit

Kali Linux Xfce(X11)  
┌┌──(alexander㉿kali)-[~/Documents/Point_Search_JLP]
└─$ ./generate_bloom
[23:01:16] P_table generated
[23:01:16] Range Start: 64 bits
[23:01:16] Range End  : 65 bits
[23:01:16] Block Width: 2^31
[23:01:16] Search Pub : 03c7290b09537769fe749210c75bc592e127204ad8b677c79f3f0c51e73b2f6d8c
[23:01:16] Settings written to file
[23:01:16] Creating bloomfilter images
[23:29:44] Writing bloom2 image to bloom2.bf
[23:29:59] Writing bloom1 image to bloom1.bf
[23:32:48] Elapsed time: (0)hours (31)minutes (32)seconds
                                                                                                                  
┌──(alexander㉿kali)-[~/Documents/Point_Search_JLP]
└─$ ./point_search  
[00:01:24] S_table generated
[00:01:24] Range Start: 64 bits
[00:01:24] Range End  : 65 bits
[00:01:24] Block Width: 2^31
[00:01:24] Search Pub : 03c7290b09537769fe749210c75bc592e127204ad8b677c79f3f0c51e73b2f6d8c
[00:01:24] Loading Bloomfilter images
[00:03:33] Search in progress...
[00:04:03] BloomFilter Hit bloom1.bf (Even Point) [Lower Range Half]
[00:04:03] False Positive
[00:05:01] BloomFilter Hit bloom1.bf (Even Point) [Lower Range Half]
[00:05:01] Private key: 24624325287185305508
[00:05:01] Elapsed time: (0)hours (1)minutes (27)seconds
</pre>
