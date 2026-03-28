import random
import secp256k1

N = 115792089237316195423570985008687907852837564279074904382605163141518161494337
a_power = 53
b_power = 54
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
