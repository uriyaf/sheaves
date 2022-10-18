from random import randrange, randint



def PrimeField(p):
    class Fp(object):
        """
        Implements the Galois field with p elements. p is required to be prime.
        All members of the field are instanciated and are unique.
        No multiplication and inverse tables are stored.
        """
        _elements = None

        def __new__(cls, value):
            #if cls._elements is None:
            #    cls._elements = [object.__new__(cls) for i in range(p)]
            #    for i in range(p):
            #        cls._elements[i].a =i 
            return cls._elements[value % p]

        def __str__(self):
            return str(self.a)

        def __repr__(self):
            return str(self)

        def __add__(self,other):
            return Fp(self.a + other.a)

        def __sub__(self,other):
            return Fp(self.a - other.a)

        def __mul__(self,other):
            return Fp(self.a * other.a)

        def __truediv__(self,other):
            return (self * other.inverse())

        def __neg__(self):
            return Fp(-self.a)

        def __ne__(self,other):
            return (self.a!=other.a)

        def __eq__(self,other):
            return (self.a==other.a)

        def inverse(self):
            return Fp(pow(self.a,p-2,p))

        def rand():
            return Fp._elements[randrange(0,p)]
            
        def enum():
            for i in range(p):
                yield Fp._elements[i]

        def one():
            return Fp._elements[1]

        def zero():
            return Fp._elements[0]
    # manually instanciate elements
    Fp._elements = [object.__new__(Fp) for i in range(p)]
    for i in range(p):
        Fp._elements[i].a =i 
    return Fp

class F2(object):
    """
    Implements the Galois field with 2 elements.
    Elements are unique.
    """
    _elements = None
    
    def __new__(cls, value):
            #if cls._elements is None:
            #    cls._elements = [object.__new__(cls) for i in range(p)]
            #    for i in range(p):
            #        cls._elements[i].a =i 
            return cls._elements[value % 2]

    def __str__(self):
        return str(self.a)

    def __repr__(self):
        return str(self)

    def __add__(self,other):
        return F2(self.a ^ other.a)

    def __sub__(self,other):
        return F2(self.a ^ other.a)

    def __mul__(self,other):
        return F2(self.a & other.a)

    def __truediv__(self,other):
        return (self * other.inverse())

    def __neg__(self):
        return F2(self.a)

    def __ne__(self,other):
        return (self.a!=other.a)

    def __eq__(self,other):
        return (self.a==other.a)

    def inverse(self):
        return F2(self.a)

    def rand():
        return F2(randrange(0,2))

    def enum():
        yield F2(0)
        yield F2(1)

    def one():
        return F2(1)

    def zero():
        return F2(0)

# manually instanciate F2

F2._elements = [object.__new__(F2), object.__new__(F2)]
F2._elements[0].a = 0
F2._elements[1].a = 1


class F4(object):
    """
    Implements the Galois field with 4 elements: F2[u | u^2 + u = 1]
    """
    def __init__(self, a, b=0):
        "a+b*u, where u^2 + u = 1"
        self.a = a
        self.b = b

    def __str__(self):
        return "F4(%d,%d)" % (self.a, self.b)

    def __repr__(self):
        return str(self)

    def __add__(self,other):
        return F4(self.a ^ other.a, self.b ^ other.b)

    def __sub__(self,other):
        return F4(self.a ^ other.a, self.b ^ other.b)

    def __mul__(self,other):
        return F4((self.a * other.a)^(self.b * other.b), (self.a * other.b) ^ (self.b * other.b) ^ (self.b * other.a))

    def __truediv__(self,other):
        return (self * other.inverse())

    def __neg__(self):
        return F4(self.a,self.b)

    def __ne__(self,other):
        return (self.a!=other.a)or(self.b!=other.b)

    def __eq__(self,other):
        return (self.a==other.a)and(self.b==other.b)

    def inverse(self):
        if self.b == 0:
            return F4(self.a,0)
        return F4(1-self.a,1)

    def rand():
        return F4(randrange(0,2),randrange(0,2))

    def enum():
        yield F4(0,0)
        yield F4(1,0)
        yield F4(0,1)
        yield F4(1,1)

    def one():
        return F4(1,0)

    def zero():
        return F4(0,0)

def bits_conv(x,len_x,y,len_y):
    res = 0
    for i in range(len_y):
        if (y>>i) & 1:
            res ^= x << i
    return res

def bits_mod(x,len_x,modulus):
    if modulus == 0:
        raise "Modulus is 0"
    len_mod = 0
    while modulus >> len_mod:
        len_mod += 1
    while len_x >= len_mod:
        if (x >> (len_x-1)) & 1:
            x ^= modulus << (len_x-len_mod)
        len_x -= 1
    return x

def F16_mult_table(modulus_pol):
    res = [[] for i in range(16)]
    for i in range(16):
        for j in range(16):
            res[i].append(bits_mod(bits_conv(i,4,j,4),7,modulus_pol))
    return res

def F16_inv_table(modulus_pol):
    prod_table = F16_mult_table(modulus_pol)
    return [0] + [prod_table[i].index(1) for i in range(1,16)]

_F16_mult_table_mod_19 = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], \
                         [0, 2, 4, 6, 8, 10, 12, 14, 3, 1, 7, 5, 11, 9, 15, 13], \
                         [0, 3, 6, 5, 12, 15, 10, 9, 11, 8, 13, 14, 7, 4, 1, 2], \
                         [0, 4, 8, 12, 3, 7, 11, 15, 6, 2, 14, 10, 5, 1, 13, 9], \
                         [0, 5, 10, 15, 7, 2, 13, 8, 14, 11, 4, 1, 9, 12, 3, 6], \
                         [0, 6, 12, 10, 11, 13, 7, 1, 5, 3, 9, 15, 14, 8, 2, 4], \
                         [0, 7, 14, 9, 15, 8, 1, 6, 13, 10, 3, 4, 2, 5, 12, 11], \
                         [0, 8, 3, 11, 6, 14, 5, 13, 12, 4, 15, 7, 10, 2, 9, 1], \
                         [0, 9, 1, 8, 2, 11, 3, 10, 4, 13, 5, 12, 6, 15, 7, 14], \
                         [0, 10, 7, 13, 14, 4, 9, 3, 15, 5, 8, 2, 1, 11, 6, 12], \
                         [0, 11, 5, 14, 10, 1, 15, 4, 7, 12, 2, 9, 13, 6, 8, 3], \
                         [0, 12, 11, 7, 5, 9, 14, 2, 10, 6, 1, 13, 15, 3, 4, 8], \
                         [0, 13, 9, 4, 1, 12, 8, 5, 2, 15, 11, 6, 3, 14, 10, 7], \
                         [0, 14, 15, 1, 13, 3, 2, 12, 9, 7, 6, 8, 4, 10, 11, 5], \
                         [0, 15, 13, 2, 9, 6, 4, 11, 1, 14, 12, 3, 8, 7, 5, 10]]

_F16_inv_table = [0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8]

class F16(object):
    """
    Implements the Galois field with 16 elements: F2[u | u^4 + u + 1 = 0].
    To change the modulus polynomial, update _modulus, _product_table and _inverse_table
    """
    _modulus = 19
    _product_table = _F16_mult_table_mod_19
    _inverse_table = _F16_inv_table
    
    def __init__(self, bits=0):
        self._bits = bits

    def __str__(self):
        return str((self._bits>>3)&1) + \
               str((self._bits>>2)&1) + \
               str((self._bits>>1)&1) + \
               str((self._bits)&1)

    def __repr__(self):
        return str(self)

    def __add__(self,other):
        return F16(self._bits ^ other._bits)

    def __sub__(self,other):
        return F16(self._bits ^ other._bits)

    def __mul__(self,other):
        return F16(F16._product_table[self._bits][other._bits])

    def __truediv__(self,other):
        return (self * other.inverse())

    def __neg__(self):
        return F16(self._bits)

    def __ne__(self,other):
        return (self._bits != other._bits)

    def __eq__(self,other):
        return (self._bits == other._bits)

    def inverse(self):
        return F16(F16._inverse_table[self._bits])

    def rand():
        return F16(randint(0,15))

    def enum():
        for i in range(16):
            yield F16(i)

    def one():
        return F16(1)

    def zero():
        return F16(0)

    def to_int(self):
        return self._bits



