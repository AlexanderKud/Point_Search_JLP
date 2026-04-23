import random
import secp256k1
import multiprocessing as mp

a_power = 49
b_power = 50
block_width = (a_power // 2) if (a_power // 2) < 32 else 32
a = 2**a_power
b = 2**b_power
print(b)
print(a + a - 2**(a_power - 1))
print(a)
print(hex(b)[2:])
print(hex(a)[2:])
pk = random.randrange(a, b)
P = secp256k1.scalar_multiplication(pk)
cpub = secp256k1.point_to_cpub(P)
print(cpub)
print(pk)
print(hex(pk))
print(secp256k1.publickey_to_address(0, True, P))
print(secp256k1.publickey_to_hash160(0, True, P))

f = open("settings.txt", "w")
f.write(f"{str(a_power)}\n")
f.write(f"{str(b_power)}\n")
f.write(f"{str(block_width)}\n")
f.write(f"{cpub}\n")
f.close()

print()

num_threads = mp.cpu_count()

partition_Size = a // num_threads
center_Num = partition_Size // 2

range_Start = a
range_Nums = []
for i in range(num_threads + 1):
    range_Nums.append(range_Start)
    range_Start += partition_Size


center_Int = 0
center_Nums = []
for i in range(num_threads):
    center_Int = range_Nums[i] + center_Num
    center_Nums.append(center_Int)

for i in range(num_threads):
    print(f'{range_Nums[i]}-{center_Nums[i]}-{range_Nums[i + 1]}')
