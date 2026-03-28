import ctypes
import platform

if platform.system() == 'Linux':
    secp256k1 = ctypes.CDLL("./secp256k1_lib.so")
elif platform.system() == 'Windows':
    secp256k1 = ctypes.CDLL("./secp256k1_lib.dll")

secp256k1.check.argtypes = None
secp256k1.check.restype = None

secp256k1.scalar_multiplication.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.scalar_multiplication.restype = None

secp256k1.point_multiplication.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
secp256k1.point_multiplication.restype = None

secp256k1.point_division.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
secp256k1.point_division.restype = None

secp256k1.double_point.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.double_point.restype = None

secp256k1.negate_point.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.negate_point.restype = None

secp256k1.add_points.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
secp256k1.add_points.restype = None

secp256k1.add_point_scalar.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
secp256k1.add_point_scalar.restype = None

secp256k1.subtract_points.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
secp256k1.subtract_points.restype = None

secp256k1.subtract_point_scalar.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
secp256k1.subtract_point_scalar.restype = None

secp256k1.increment_point.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.increment_point.restype = None

secp256k1.decrement_point.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.decrement_point.restype = None

secp256k1.point_on_curve.argtypes = [ctypes.c_char_p]
secp256k1.point_on_curve.restype = ctypes.c_bool

secp256k1.get_y.argtypes = [ctypes.c_char_p, ctypes.c_bool, ctypes.c_char_p]
secp256k1.get_y.restype = None

secp256k1.privatekey_to_hash160.argtypes = [ctypes.c_int, ctypes.c_bool, ctypes.c_char_p, ctypes.c_char_p]
secp256k1.privatekey_to_hash160.restype = None

secp256k1.publickey_to_hash160.argtypes = [ctypes.c_int, ctypes.c_bool, ctypes.c_char_p, ctypes.c_char_p]
secp256k1.publickey_to_hash160.restype = None

secp256k1.privatekey_to_uwif.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.privatekey_to_uwif.restype = None

secp256k1.privatekey_to_cwif.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.privatekey_to_cwif.restype = None

secp256k1.wif_to_privatekey.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.wif_to_privatekey.restype = None

secp256k1.privatekey_to_address.argtypes = [ctypes.c_int, ctypes.c_bool, ctypes.c_char_p, ctypes.c_char_p]
secp256k1.privatekey_to_address.restype = None

secp256k1.publickey_to_address.argtypes = [ctypes.c_int, ctypes.c_bool, ctypes.c_char_p, ctypes.c_char_p]
secp256k1.publickey_to_address.restype = None

secp256k1.privatekey_to_bech32_address.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.privatekey_to_bech32_address.restype = None

secp256k1.publickey_to_bech32_address.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.publickey_to_bech32_address.restype = None

secp256k1.publickey_to_bech32_p2wsh_address.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.publickey_to_bech32_p2wsh_address.restype = None

secp256k1.hash160_to_address.argtypes = [ctypes.c_int, ctypes.c_bool, ctypes.c_char_p, ctypes.c_char_p]
secp256k1.hash160_to_address.restype = None

secp256k1.publickey_to_point.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.publickey_to_point.restype = None

secp256k1.p2pkh_address_to_hash160.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
secp256k1.p2pkh_address_to_hash160.restype = None

secp256k1.init_bloom.argtypes = [ctypes.c_int, ctypes.c_ulonglong, ctypes.c_double]
secp256k1.init_bloom.restype = None

secp256k1.bloom_info.argtypes = [ctypes.c_int]
secp256k1.bloom_info.restype = None

secp256k1.bloom_save.argtypes = [ctypes.c_int, ctypes.c_char_p]
secp256k1.bloom_save.restype = None

secp256k1.bloom_load.argtypes = [ctypes.c_int, ctypes.c_char_p]
secp256k1.bloom_load.restype = None

secp256k1.bloom_add.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int]
secp256k1.bloom_add.restype = None

secp256k1.bloom_check.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int]
secp256k1.bloom_check.restype = ctypes.c_int

N       = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
lambda1 = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
lambda2 = 0xac9c52b33fa3cf1f5ad9e3fd77ed9ba4a880b9fc8ec739c2e0cfc810b51283ce
P       = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
beta1   = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
beta2   = 0x851695d49a83f8ef919bb86153cbcb16630fb68aed0a766a3ec693d68e6afa40
    
secp256k1.Init()

def check():
    secp256k1.check()
    
def scalar_multiplication(pk):
    pvk = (pk % N).to_bytes(32, 'big')
    res = bytes(65)
    secp256k1.scalar_multiplication(pvk, res)
    return res

def point_multiplication(p, pk):
    pvk = (pk % N).to_bytes(32, 'big')
    res = bytes(65)
    secp256k1.point_multiplication(p, pvk, res)
    return res
   
def point_division(p, pk):
    pvk = (pk % N).to_bytes(32, 'big')
    res = bytes(65)
    secp256k1.point_division(p, pvk, res)
    return res

def point_to_upub(pBytes):
    return pBytes.hex()

def point_to_cpub(pBytes):
    ph = pBytes.hex()
    last_byte = int(ph[128:130], 16)
    prefix = '02' if (last_byte & 1) == 0 else '03'
    return prefix + ph[2:66]

def double_point(pBytes):
    res = bytes(65)
    secp256k1.double_point(pBytes, res)
    return res

def negate_point(pBytes):
    res = bytes(65)
    secp256k1.negate_point(pBytes, res)
    return res

def point_endo1(pBytes):
    x_int = (int.from_bytes(pBytes[1:33], 'big') * beta1) % P
    x_bytes = x_int.to_bytes(32, 'big')
    return pBytes[:1] + x_bytes + pBytes[33:]

def point_endo2(pBytes):
    x_int = (int.from_bytes(pBytes[1:33], 'big') * beta2) % P
    x_bytes = x_int.to_bytes(32, 'big')
    return pBytes[:1] + x_bytes + pBytes[33:]
    
def priv_endo1(pk):
    return (pk * lambda1) % N

def priv_endo2(pk):
    return (pk * lambda2) % N
    
def add_points(p1, p2):
    res = bytes(65)
    secp256k1.add_points(p1, p2, res)
    return res

def add_point_scalar(p, pk):
    pvk = (pk % N).to_bytes(32, 'big')
    res = bytes(65)
    secp256k1.add_point_scalar(p, pvk, res)
    return res

def subtract_points(p1, p2):
    res = bytes(65)
    secp256k1.subtract_points(p1, p2, res)
    return res

def subtract_point_scalar(p, pk):
    pvk = (pk % N).to_bytes(32, 'big')
    res = bytes(65)
    secp256k1.subtract_point_scalar(p, pvk, res)
    return res

def increment_point(pBytes):
    res = bytes(65)
    secp256k1.increment_point(pBytes, res)
    return res

def decrement_point(pBytes):
    res = bytes(65)
    secp256k1.decrement_point(pBytes, res)
    return res

def point_on_curve(pBytes):
    return secp256k1.point_on_curve(pBytes)
    
# addr_type = 0 [p2pkh] compressed=False(uncompressed) compressed=True(compressed),  1 [p2sh] compressed=True
def privatekey_to_hash160(addr_type, compressed, pk):
    pvk = (pk % N).to_bytes(32, 'big')
    res = bytes(20)
    secp256k1.privatekey_to_hash160(addr_type, compressed, pvk, res)
    return res.hex()
    
# addr_type = 0 [p2pkh] compressed=False(uncompressed) compressed=True(compressed),  1 [p2sh] compressed=True    
def publickey_to_hash160(addr_type, compressed, pBytes):
    res = bytes(20)
    secp256k1.publickey_to_hash160(addr_type, compressed, pBytes, res)
    return res.hex()

def privatekey_to_uwif(pk):
    pvk = (pk % N).to_bytes(32, 'big')
    res = bytes(51)
    secp256k1.privatekey_to_uwif(pvk, res)
    return res.decode('utf-8')

def privatekey_to_cwif(pk):
    pvk = (pk % N).to_bytes(32, 'big')
    res = bytes(52)
    secp256k1.privatekey_to_cwif(pvk, res)
    return res.decode('utf-8')

def wif_to_privatekey(wif):
    res = bytes(32)
    secp256k1.wif_to_privatekey(wif.encode('utf-8'), res)
    return int.from_bytes(res, 'big')

# addr_type = 0 [p2pkh] compressed=False(uncompressed) compressed=True(compressed),  1 [p2sh] compressed=True
def privatekey_to_address(addr_type, compressed, pk):
    pvk = (pk % N).to_bytes(32, 'big')
    res = bytes(35)
    secp256k1.privatekey_to_address(addr_type, compressed, pvk, res)
    return res.rstrip(b'\x00').decode('utf-8')

# addr_type = 0 [p2pkh] compressed=False(uncompressed) compressed=True(compressed),  1 [p2sh] compressed=True    
def publickey_to_address(addr_type, compressed, p):
    res = bytes(35)
    secp256k1.publickey_to_address(addr_type, compressed, p, res)
    return res.rstrip(b'\x00').decode('utf-8')

def privatekey_to_bech32_address(pk):
    pvk = (pk % N).to_bytes(32, 'big')
    res = bytes(42)
    secp256k1.privatekey_to_bech32_address(pvk, res)
    return res.decode('utf-8')

def publickey_to_bech32_address(p):
    res = bytes(42)
    secp256k1.publickey_to_bech32_address(p, res)
    return res.decode('utf-8')

def publickey_to_bech32_p2wsh_address(p):
    res = bytes(62)
    secp256k1.publickey_to_bech32_p2wsh_address(p, res)
    return res.decode('utf-8')

def hash160_to_address(addr_type, compressed, hash160):
    res = bytes(42)
    secp256k1.hash160_to_address(addr_type, compressed, bytes.fromhex(hash160), res)
    return res.rstrip(b'\x00').decode('utf-8')

def publickey_to_point(pub):
    x = pub[2:66]
    if len(pub) == 66:
        res = bytes(32)
        secp256k1.get_y(bytes.fromhex(x), int(pub[:2], 16) % 2 == 0, res)
        y = res.hex()
    else:
        y = pub[66:]
    return bytes.fromhex('04' + x + y)

def p2pkh_address_to_hash160(address):
    res = bytes(25)
    secp256k1.p2pkh_address_to_hash160(address.encode('utf-8'), res)
    return res.hex()[2:42]
    
def init_bloom(index, entries, error):
    secp256k1.init_bloom(index, entries, error)

def bloom_info(index):
    secp256k1.bloom_info(index)

def bloom_save(index, filename):
    secp256k1.bloom_save(index, filename.encode('utf-8'))
    
def bloom_load(index, filename):
    secp256k1.bloom_load(index, filename.encode('utf-8'))

def bloom_add(index, item):
    if type(item) != bytes: item = str(item).encode('utf-8')
    secp256k1.bloom_add(index, item, len(item))

def bloom_check(index, item):
    if type(item) != bytes: item = str(item).encode('utf-8')
    return secp256k1.bloom_check(index, item, len(item))

def bloom_add_bytes(index, item):
    secp256k1.bloom_add(index, item, len(item))

def bloom_check_bytes(index, item):
    return secp256k1.bloom_check(index, item, len(item))
