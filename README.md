<pre>[alexander@alexander-home Point_Search_CPP]$ ./generate_bloom
[00:16:21] P_table generated
[00:16:21] Range Start: 50 bits
[00:16:21] Range End  : 51 bits 
[00:16:21] Block Width: 2^25
[00:16:21] Search Pub : 037fcccd965224a7d299ec4f14c1506b35b3583d3fdae17a69f9b6fec3dd2089a2
[00:16:21] Settings written to file
[00:16:21] Creating BloomFile1
[00:16:21] Creating BloomFile2
[00:25:31] Writing BloomFile1 to bloom1.bf
[00:25:33] Writing BloomFile2 to bloom2.bf
[00:25:34] Elapsed time: 552.742s

[alexander@alexander-home Point_Search_CPP]$ ./point_search
[00:25:39] S_table generated
[00:25:39] Range Start: 50 bits
[00:25:39] Range End  : 51 bits
[00:25:39] Block Width: 2^25
[00:25:39] Search Pub : 037fcccd965224a7d299ec4f14c1506b35b3583d3fdae17a69f9b6fec3dd2089a2
[00:25:39] Loading Bloomfilter bloom1.bf
[00:25:39] Loading Bloomfilter bloom2.bf
[00:25:39] Search in progress...
[00:26:45] BloomFilter Hit bloom1.bf (Even Point) [Lower Range Half]
[00:26:45] Privatekey: 1295403161015294
[00:26:45] Elapsed time: 66.0754s
</pre>
