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
- calculating just x coordinate for the batch
- calculating x,y for the last of the batch entry (used as the next startPoint)

point_search.cpp
- batch addition
- batch inversion
- calculating x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)
- bloom check by x coordinate uint64_t bits64[3] part, computing y coordinate only if there is a hit

Kali Linux Xfce(X11)  
┌──(alexander㉿kali)-[~/Documents/Point_Search_JLP_Test]
└─$ ./generate_bloom
[00:53:53] P_table generated
[00:53:53] Range Start: 59 bits
[00:53:53] Range End  : 60 bits
[00:53:53] Block Width: 2^30
[00:53:53] Search Pub : 0386d42f693530d42401a660259c79a74796db05e1ebe8bef5727c535a1f45df80
[00:53:53] Settings written to file
[00:53:53] Creating bloomfilter images
[01:07:49] Writing bloom2 image to bloom2.bf
[01:07:57] Writing bloom1 image to bloom1.bf
[01:08:58] Elapsed time: (0)hours (15)minutes (6)seconds
                                                                                                                  
┌──(alexander㉿kali)-[~/Documents/Point_Search_JLP_Test]
└─$ ./point_search  
[01:09:09] S_table generated
[01:09:09] Range Start: 59 bits
[01:09:09] Range End  : 60 bits
[01:09:09] Block Width: 2^30
[01:09:09] Search Pub : 0386d42f693530d42401a660259c79a74796db05e1ebe8bef5727c535a1f45df80
[01:09:09] Loading Bloomfilter images
[01:09:16] Search in progress...
[01:09:22] BloomFilter Hit bloom2.bf (Odd Point) [Higher Range Half]
[01:09:22] Private key: 910788673462129357
[01:09:22] Elapsed time: (0)hours (0)minutes (6)seconds


</pre>
