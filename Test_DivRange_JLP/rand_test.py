import random
from bit import Key, format

N = 115792089237316195423570985008687907852837564279074904382605163141518161494337
a_power = 39
b_power = 40
a = 2**a_power
b = 2**b_power
print(b)
print(a + a - 2**(a_power - 1))
print(a)
print(hex(b)[2:])
print(hex(a)[2:])
pk = random.randrange(a, b)
keyC = Key.from_int(pk)
print(keyC.public_key.hex())
print(pk)
print(hex(pk))
print(keyC.address)
print((format.address_to_public_key_hash(keyC.address)).hex())
