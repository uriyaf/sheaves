
from itertools import product

class Z2n(object):
    """
    Implements the Galois field with p elements. p is required to be prime.
    All members of the field are instanciated and are unique.
    No multiplication and inverse tables are stored.
    """

    def __init__(self,value, bits=1):
        """
        modulus must be a power of 2.
        """
        self.a = value & ((1 << bits) - 1)
        self._bits = bits

    def __str__(self):
        return str(self.a)

    def __repr__(self):
        return str(self)

    def __add__(self,other):
        return Z2n(self.a + other.a, self._bits)

    def __sub__(self,other):
        return Z2n(self.a - other.a, self._bits)

    def __mul__(self,other):
        return Z2n(self.a * other.a, self._bits)

    def __truediv__(self,other):
        return self * other.inverse()

    def __neg__(self):
        return Z2n(-self.a, self._bits)

    def __ne__(self,other):
        return (self.a!=other.a)

    def __eq__(self,other):
        return (self.a==other.a)

    def inverse(self):
        return Z2n(pow(self.a, 1 << (self._bits-2) - 1, 1 << self._bits), self._bits)

    def get_lifts(self):
        b = self._bits
        a = self.a
        yield Z2n(self.a, b+1)
        yield Z2n(self.a + (1<<b), b+1)


class MatrixZ2n(object):
    """
    A matrix over Z2n.
    """

    def __init__(self, rows, cols, mat=None, bits=None):
        """
        @rows is the number of rows,
        @cols is the number of columns,
        @mat is a list of lists representing the matrix.
        If @mat is omitted, a zero matrix is created.
        """
        self._cols = cols
        self._rows = rows
        if mat is None:
            assert bits is not None
            zero = Z2n(0,bits)
            self._mat = [[zero for i in range(cols)] for j in range(rows)]
        else:
            self._mat = [r[:] for r in mat]
        if bits is None:
            bits = mat[0][0]._bits
        self._bits = bits

    def rows(self):
        return self._rows

    def cols(self):
        return self._cols

    def get_entry(self,i,j):
        """
        Returns the (i,j) entry.
        """
        return self._mat[i][j]

    def __mul__(self, other):
        res = MatrixZ2n(self._rows, other._cols, mat=None, bits=self._bits)
        for i in range(self._rows):
            for j in range(other._cols):
                for k in range(self._cols):
                    res._mat[i][j] = res._mat[i][j] + self._mat[i][k] * other._mat[k][j]
        return res

    def __repr__(self):
        rows = self._rows
        cols = self._cols 
        res = "%d-by-%d matrix:\n" % (rows, cols)
        for i in range(rows):
            res += "[" + (",".join(str(self._mat[i][j]) for j in range(cols))) + "]"
            if i < rows:
                res += "\n"
        return res

    def get_lifts(self):
        b = self._bits
        L = [0, 1<<b]
        for m in product(L, repeat=self._rows * self._cols):
            mat = [[Z2n(x.a + m[i*self._rows + j],b+1) for j,x in enumerate(r)] for i,r in enumerate(self._mat)]
            yield MatrixZ2n(self._rows, self._cols, mat, b+1)

    def __eq__(self, other):
        if self._rows != other._rows:
            return False
        if self._cols != other._cols:
            return False
        for i in range(self.rows()):
            for j in range(self.cols()):
                if self._mat[i][j] != other._mat[i][j]:
                    return False
        return True

    def __neq__(self, other):
        if self._rows != other._rows:
            return True
        if self._cols != other._cols:
            return True
        for i in range(self.rows()):
            for j in range(self.cols()):
                if self._mat[i][j] != other._mat[i][j]:
                    return True
        return False

def solve_eq(max_stage=3, begin_stage=1, begin_sols=None, max_sols=30000):
    sols = begin_sols
    stage = begin_stage
    if sols is None:
        I = MatrixZ2n(3, 3, [[Z2n(1,1), Z2n(0,1), Z2n(0,1)], \
                         [Z2n(0,1), Z2n(1,1), Z2n(0,1)],
                         [Z2n(0,1), Z2n(0,1), Z2n(1,1)]])
        sols = [I]
        stage = 1
    X = None
    while stage < max_stage:
        stage += 1
        #X = MatrixZ2n(3, 3, [[Z2n(1,stage), Z2n(0,stage), Z2n(0,stage)], \
        #                     [Z2n(0,stage), Z2n(-1,stage), Z2n(0,stage)],
        #                     [Z2n(0,stage), Z2n(0,stage), Z2n(-1,stage)]])
        X = MatrixZ2n(3, 3, [[Z2n(-1,stage), Z2n(0,stage), Z2n(0,stage)], \
                             [Z2n(0,stage), Z2n(0,stage), Z2n(-1,stage)],
                             [Z2n(0,stage), Z2n(-1,stage), Z2n(0,stage)]])
        print("stage:", stage)
        new_sols = []
        for S in sols:
            for L in S.get_lifts():
                L2 = L*L
                L3 = L2*L
                if L3*X == X*L3*L2: # X == L*X*L*L2*L2:
                    new_sols.append(L)
            if len(new_sols) > max_sols:
                break
        sols = new_sols
        print("solutions:", len(sols))
    return sols, X
    
def solve_eq_adhoc(max_stage=3, begin_stage=1, begin_sols=None, max_sols=30000):
    sols = begin_sols
    stage = begin_stage
    if sols is None:
        I = MatrixZ2n(3, 3, [[Z2n(1,1), Z2n(0,1), Z2n(0,1)], \
                         [Z2n(0,1), Z2n(1,1), Z2n(0,1)],
                         [Z2n(0,1), Z2n(0,1), Z2n(1,1)]])
        sols = [I]
        stage = 1
    X = None
    while stage < max_stage:
        stage += 1
        #X = MatrixZ2n(3, 3, [[Z2n(1,stage), Z2n(0,stage), Z2n(0,stage)], \
        #                     [Z2n(0,stage), Z2n(-1,stage), Z2n(0,stage)],
        #                     [Z2n(0,stage), Z2n(0,stage), Z2n(-1,stage)]])
        X = MatrixZ2n(3, 3, [[Z2n(-1,stage), Z2n(0,stage), Z2n(0,stage)], \
                             [Z2n(0,stage), Z2n(0,stage), Z2n(-1,stage)],
                             [Z2n(0,stage), Z2n(-1,stage), Z2n(0,stage)]])
        print("stage:", stage)
        new_sols = []
        for S in sols:
            for L in S.get_lifts():
                L2 = L*L
                if X == L*X*L*L2*L2:
                    if L*X*L2*X*L != X*X:
                        new_sols.append(L)
                        print(L)
            if len(new_sols) > max_sols:
                break
        sols = new_sols
        print("solutions:", len(sols))
    return sols, X
