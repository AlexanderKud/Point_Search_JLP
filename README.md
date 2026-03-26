<pre>
generate_bloom.cpp
- batch addition
- batch inversion
- calculating x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)

point_search.cpp
- batch addition
- batch inversion
- calculating x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)

[alexander@alexander-manjaro Point_Search_JLP]$ ./generate_bloom
[18:59:05] P_table generated
[18:59:05] Range Start: 65 bits
[18:59:05] Range End  : 66 bits
[18:59:05] Block Width: 2^31
[18:59:05] Search Pub : 0217d9823862bf7a648e3d88f9685aec65be263b15a56424078ddce4ca006e012e
[18:59:05] Settings written to file
[18:59:05] Bloomfilter Size : 8.00 GB Total(2 blooms): 16.00 GB
[18:59:05] Creating bloomfilter images
[19:06:10] Writing BloomFilter to bloom1.bf
[19:06:11] Writing BloomFilter to bloom2.bf
[19:08:12] Elapsed time: (0)hours (9)minutes (7)seconds

[alexander@alexander-manjaro Point_Search_JLP]$ ./point_search
[19:10:02] S_table generated
[19:10:02] Range Start: 65 bits
[19:10:02] Range End  : 66 bits
[19:10:02] Block Width: 2^31
[19:10:02] Search Pub : 0217d9823862bf7a648e3d88f9685aec65be263b15a56424078ddce4ca006e012e
[19:10:02] Loading Bloomfilter images
[19:10:20] Search in progress...
[19:14:46] Save Data written to settings1.txt
[19:14:47] Save Data written to settings2.txt
[19:18:04] BloomFilter Hit bloom1.bf (Even Point) [Lower Range Half]
[19:18:04] Private key: 44179411046104651556
[19:18:04] Elapsed time: (0)hours (7)minutes (43)seconds

Point_Search2_JLP
[alexander@alexander-manjaro Point_Search2_JLP]$ ./bloom
[13:33:26] P_table generated
[13:33:26] Range Start: 69 bits
[13:33:26] Range End  : 70 bits
[13:33:26] Block Width: 2^32
[13:33:26] Search Pub : 03e2929572e326c95ab0a6da7324d2410dd6328f27040a876fd880a399aab8da3a
[13:33:26] Settings written to file
[13:33:26] Bloomfilter Size : 16.00 GB
[13:33:26] Creating bloomfilter image
[13:46:37] Writing BloomFilter to bloom.bf
[13:48:57] Elapsed time: (0)hours (15)minutes (31)seconds

[alexander@alexander-manjaro Point_Search2_JLP]$ ./search
[13:49:16] S_table generated
[13:49:16] Range Start: 69 bits
[13:49:16] Range End  : 70 bits
[13:49:16] Block Width: 2^32
[13:49:16] Search Pub : 03e2929572e326c95ab0a6da7324d2410dd6328f27040a876fd880a399aab8da3a
[13:49:16] Loading Bloomfilter image
[13:49:39] Search in progress...
[13:53:31] Save Data written to stride_sum.txt
[13:57:21] Save Data written to stride_sum.txt
[14:01:10] Save Data written to stride_sum.txt
[14:04:59] Save Data written to stride_sum.txt
[14:08:55] Save Data written to stride_sum.txt
[14:12:51] Save Data written to stride_sum.txt
[14:13:54] Found by Thread -> 1
[14:13:54] Private key: 700782406099173415617
[14:13:54] Elapsed time: (0)hours (24)minutes (15)seconds
</pre>
