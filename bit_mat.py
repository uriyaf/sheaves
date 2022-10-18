
from base_mat import *
from random import randint
from itertools import chain

class BitVec(BaseVector):
    """
    Vector of bits.
    """

    def __init__(self, size, int_vec = None, bit_vec = None, block_size=3):
        """
        @size is the length of the vector.
        Specify @int_vec to initialize from an int list.
        Specify @bits (of the given block size) to intialize form list representing a bit vector.
        """
        self._len = size
        self._bs = block_size
        if bit_vec is None:
            self._vec = [0 for i in range((size+block_size-1) // block_size)]
            if int_vec is not None:
                for i,v in enumerate(int_vec):
                    self.set_entry(i,v)
        else:
            self._vec = bit_vec

    # implementation of methods of BaseVector

    def get_entry(self,i):
        """
        Returns the i-th entry of the vector.
        """
        return (self._vec[i // self._bs] >> (i % self._bs)) & 1

    def set_entry(self,i,v):
        """
        Returns the i-th entry to v
        """
        v = v & 1
        of = i % self._bs
        b = i // self._bs
        self._vec[b] |= (1 << of)
        self._vec[b] ^= ((1-v) << of)

    def vec_class(self):
        return BitVecClass(self._bs)

    def __len__(self):
        return self._len

    def clone(self):
        return BitVec(self._len, bit_vec = self._vec[:], block_size = self._bs)

    def scale(self, scalar):
        """
        Returns the vector obtained from scaling @self by the given scalar.
        """
        scalar = scalar & 1
        if scalar == 1:
            return self.clone()
        else:
            return BitVec(self._len, block_size = self._bs)

    # new methods:

    def rand(size, block_size=120):
        b = (size+block_size-1)//block_size
        bit_vec = [randint(0, (1 << block_size)-1) for i in range(size//block_size)]
        of = size % block_size
        if of > 0:
            bit_vec.append(randint(0, (1 << of) - 1))
        return BitVec(size, bit_vec = bit_vec, block_size = block_size)

    # overriding inherited methods:

    def __add__(self, other):
        if self._bs != other._bs:
            return BaseVector.__add__(self, other)
        bit_vec = [x^y for x,y in zip(self._vec, other._vec)]
        return BitVec(self._len, bit_vec=bit_vec, block_size=self._bs)

    def __sub__(self, other):
        if self._bs != other._bs:
            return BaseVector.__sub__(self, other)  
        bit_vec = [x^y for x,y in zip(self._vec, other._vec)]
        return BitVec(self._len, bit_vec = bit_vec, block_size = self._bs)

    def __neg__(self):
        return BitVec(self._len, bit_vec = self._vec[:], block_size = self._bs)

    def __eq__(self,other):
        if self._bs != other._bs:
            return BaseVector.__eq__(self, other)
        if len(self) != len(other):
            return False
        return self._vec == other._vec

    def __neq__(self,other):
        if self._bs != other._bs:
            return BaseVector.__neq__(self, other)
        if len(self) != len(other):
            return True
        return self._vec != other._vec

    def is_zero(self):
        for v in self._vec:
            if v != 0:
                return False
        return True

    def __mul__(self,other):
        """
        Returns the dot-product of @self and @other.
        """
        if self._bs != other._bs:
            return BaseVector.__mul__(self, other) & 1
        res = 0
        for x,y in zip(self._vec,other._vec):
            res ^= x & y
        bs = self._bs
        while bs > 1:
            half = (bs + 1) // 2
            res ^= res >> half
            res &= (1 << half) - 1
            bs = half
        return res

    def __repr__(self):
        return "".join(str(self.get_entry(i)) for i in range(self._len))
    

class BitVecClass(BaseVecClass):
    """
    A metaclass for the vector type of BitVec.
    This class holds methods which are not instant specific.
    It should be a singleton,
    but I did not have time to implement a proper singleton design pattern, though.
    """

    def __init__(self, block_size):
        self._bs = block_size

    def zero_vec(self,n):
        """
        Returns a zero vector of length n.
        """
        return BitVec(n, block_size = self._bs)

    def rand(self,n):
        """
        Returns a random vector of length n.
        """
        return BitVec.rand(n, self._bs)

    def __eq__(self,other):
        return self._bs == other._bs

    # functions you get for free:

    def __neq__(self,other):
        return self._bs != other._bs






class BitMat(BaseMatrix):
    """
    Implements a bit matrix.
    """

    def __init__(self, rows, columns, mat=None, block_size=60, bit_mat=None):
        """
        Use mat to initialize from a 0s and 1s matrix.
        Use bit_mat to initialize from a bit matrix in block form (each entry stands for 'block_size' columns).
        Do not use both.
        If neither mat nor bit_mat is specified, a zero matrix is produced.
        """
        self._rows = rows
        self._cols = columns
        self._bs = block_size
        self._bcols = (columns+block_size-1) // block_size
        if bit_mat is None:
            self._mat = [[0 for i in range(self._bcols)] for j in range(self._rows)]
        else:
            self._mat = bit_mat
        if mat is not None:
            for r in range(rows):
                for c in range(columns):
                    self.set_bit(r,c,mat[r][c])

    # must implement for ConstBaseMat

    def rows(self):
        return self._rows

    def cols(self):
        return self._cols

    def get_entry(self,i,j):
        """
        Returns the (i,j) entry.
        """
        return (self._mat[i][j//self._bs] >> (j % self._bs)) & 1

    def Gauss_elim_extended_extras(self,extra):
        """
        @extra is another matrix with the same number of rows as @self.
        Performs Gauss elimination on @self.
        The row operations done on @self are also performed on a copy @extra.
        Returns the Gauss elimination of @self, a list of pairs with the leading entries, and the resulting @extra.
        """
        if self._rows != extra.rows():
            raise "Number of rows do not match."
        cols = self._cols
        rows = self._rows
        bcols = self._bcols
        if (rows == 0) or (cols == 0):
            return BitMat(0, cols, block_size = self._bs), [], extra.clone()
        E = self._bit_array_clone()
        B = extra._bit_array_clone()
        bcols_B = (extra.cols() + extra.get_block_size() - 1) // extra.get_block_size()
        leading = []
        i_b = 0  # leading column block 
        j = 0    # leading row
        while i_b < bcols:
            max_of = cols - i_b*self._bs
            i_of = 0 # leading column offset 
            while i_of < max_of:
                s = j
                while s < rows and (E[s][i_b] >> i_of) & 1 == 0 :
                    s += 1
                if s < rows:
                    leading.append((j, i_b*self._bs + i_of))
                    E[s], E[j] = E[j], E[s]
                    B[s], B[j] = B[j], B[s]
                    for l in chain(range(0,j),range(j+1,rows)):
                        if (E[l][i_b] >> i_of) & 1 == 0:
                            continue
                        for k in range(i_b,bcols):
                            E[l][k] = E[l][k] ^ E[j][k]
                        for k in range(bcols_B):
                            B[l][k] = B[l][k] ^ B[j][k]
                    j = j+1
                i_of = i_of+1
            i_b = i_b+1
        return BitMat(j,cols, block_size = self._bs, bit_mat = E[:j]), \
               leading, \
               BitMat(rows, extra.cols(), block_size = extra.get_block_size(), bit_mat = B)

    def clone(self):
        clone_mat = [r[:] for r in self._mat]
        return BitMat(self._rows, self._cols, block_size = self._bs, bit_mat = clone_mat)

    def mat_class(self):
        """
        Returns the matrix class of the given object.
        Use the matrix class to access various functions which are not instant specific, e.g.,
        constructors of that class.
        (Note: this is not the same as type().)
        """
        return BitMatClass(self._bs)

    # must implement for BaseMat

    def set_entry(self,i,j,v):
        """
        Sets entry (i,j) to v.
        """
        v = v & 1
        of = j % self._bs
        col = j // self._bs
        self._mat[i][col] |= (1 << of)
        self._mat[i][col] ^= ((1-v) << of)

    def set_to_zero(self,i,j):
        """
        Sets entry (i,j) to 0.
        """
        self._mat[i][j // self._bs] &= (1 << self._bs) - 1 - (1 << (j % self._bs))

    def set_to_one(self,i,j):
        """
        Sets entry (i,j) to 1.
        """
        self._mat[i][j // self._bs] |= (1 << (j % self._bs))

    def add_zero_rows(self, rows):
        """
        Adds @rows rows of zeroes to @self.
        """
        for i in range(rows):
            self._mat.append([0 for j in range(self._bcols)])
        self._rows += rows

    def append_rows(self, other):
        """
        Adds the rows of @other at the end of @self.
        """
        if self._cols != other.cols():
            raise "Number of columns incompatible."
        if self._bs == other.get_block_size():
            self._mat = self._mat + other._bit_array_clone()
            self._rows += other.rows()
        else:
            rows = self._rows
            other_rows = other.rows()
            self.add_zero_rows(other_rows)
            for i in range(other_rows):
                for j in range(other.cols()):
                    self.set_entry(i + rows, j, other.get_entry(i,j))

    # functions for the bit matrices family of classes

    def _bit_array_clone(self):
        """
        Returns a copy of the internal 2-dim array. Not meant for users.
        """
        return [r[:] for r in self._mat]

    def get_block_size(self):
        """
        Returns the size of each bit block.
        """
        return self._bs

    def _bit_array(self):
        """
        If mutable, returns actual inner bit array.
        If immutable, returns a copy of inner bit array.
        """
        return self._mat

    # overriding implementations of the base classes ConstBaseMatrix and BaseMatrix

    def get_row(self,i):
        """
        Returns a vector representing the i-th row.
        """
        return BitVec(self._cols, bit_vec = self._mat[i][:], block_size = self._bs)

    def get_col(self,j):
        """
        Returns a vector representing the i-th column.
        """
        b = j // self._bs
        of = j % self._bs
        col = [(self._mat[i][b] >> of) & 1 for i in range(self._rows)]
        return BitVec(self._rows, int_vec = col, block_size = self._bs)

    def get_col_as(self,j,vec_class):
        """
        Returns a vector representing the i-th column. The result has the given vector class.
        """
        b = j // self._bs
        of = j % self._bs
        rows = self._rows
        res = vec_class.zero_vec(rows)
        for i in range(rows):
                res.set_entry(i, (self._mat[i][b] >> of) & 1)
        return res

    def Gauss_elim_extended(self):
        """
        Returns the Gauss elimination of the matrix and a list of the leading entries (pairs of coordinates)
        """
        cols = self._cols
        rows = self._rows
        bcols = self._bcols
        if (rows == 0) or (cols == 0):
            return BitMat(0, cols, block_size = self._bs), []
        E = self._bit_array_clone()
        leading = []
        i_b = 0  # leading column block 
        j = 0    # leading row
        while i_b < bcols:
            max_of = cols - i_b*self._bs
            i_of = 0 # leading column offset 
            while i_of < max_of:
                s = j
                while s < rows and (E[s][i_b] >> i_of) & 1 == 0 :
                    s += 1
                if s < rows:
                    leading.append((j, i_b*self._bs + i_of))
                    E[s], E[j] = E[j], E[s]
                    for l in chain(range(0,j),range(j+1,rows)):
                        if (E[l][i_b] >> i_of) & 1 == 0:
                            continue
                        for k in range(i_b,bcols):
                            E[l][k] = E[l][k] ^ E[j][k]
                    j = j+1
                i_of = i_of+1
            i_b = i_b+1
        return BitMat(j,cols, block_size = self._bs, bit_mat = E[:j]), \
               leading
        
    def Gauss_elim(self):
        """
        Returns the Gauss elimination of the matrix
        """
        cols = self._cols
        rows = self._rows
        bcols = self._bcols
        if (rows == 0) or (cols == 0):
            return BitMat(0, cols, block_size = self._bs)
        E = self._bit_array_clone()
        i_b = 0  # leading column block 
        j = 0    # leading row
        while i_b < bcols:
            max_of = cols - i_b*self._bs
            i_of = 0 # leading column offset 
            while i_of < max_of:
                s = j
                while s < rows and (E[s][i_b] >> i_of) & 1 == 0 :
                    s += 1
                if s < rows:
                    E[s], E[j] = E[j], E[s]
                    for l in chain(range(0,j),range(j+1,rows)):
                        if (E[l][i_b] >> i_of) & 1 == 0:
                            continue
                        for k in range(i_b,bcols):
                            E[l][k] = E[l][k] ^ E[j][k]
                    j = j+1
                i_of = i_of+1
            i_b = i_b+1
        return BitMat(j, cols, block_size = self._bs, bit_mat = E[:j])

    def __add__(self,other):
        mat_class = other.mat_class()
        if isinstance(mat_class, BitMatClass):
            if self._bs == mat_class._bs:
                if isinstance(other, BitMat):
                    m1 = self._bit_array_clone()
                    m2 = other._mat
                    for i in range(self._rows):
                        for j in range(self._bcols):
                            m1[i][j] ^= m2[i][j]
                    return BitMat(self._rows, self._cols, block_size = self._bs, bit_mat = m1) 
                else:
                    m1 = self._bit_array_clone()
                    vec_class = BitVecClass(self._bs)
                    for i in range(self._rows):
                        row = other.get_row_as(i, vec_class)
                        for j in range(self._bcols):
                            m1[i][j] ^= row._vec[j]
                    return BitMat(self._rows, self._cols, block_size = self._bs, bit_mat = m1) 
            else:
                return ConstBaseMatrix.__add__(self,other)
        return ConstBaseMatrix.__add__(self,other)

    def __sub__(self,other):
        mat_class = other.mat_class()
        if isinstance(mat_class, BitMatClass):
            if self._bs == mat_class._bs:
                if isinstance(other, BitMat):
                    m1 = self._bit_array_clone()
                    m2 = other._mat
                    for i in range(self._rows):
                        for j in range(self._bcols):
                            m1[i][j] ^= m2[i][j]
                    return BitMat(self._rows, self._cols, block_size = self._bs, bit_mat = m1) 
                else:
                    m1 = self._bit_array_clone()
                    vec_class = BitVecClass(self._bs)
                    for i in range(self._rows):
                        row = other.get_row_as(i, vec_class)
                        for j in range(self._bcols):
                            m1[i][j] ^= row._vec[j]
                    return BitMat(self._rows, self._cols, block_size = self._bs, bit_mat = m1) 
            else:
                return ConstBaseMatrix.__add__(self,other)
        return ConstBaseMatrix.__add__(self,other)

    def __neg__(self):
        return self.clone()

    def __eq__(self,other):
        mat_class = other.mat_class()
        if isinstance(mat_class, BitMatClass):
            if self._bs == mat_class._bs:
                if self._rows != other.rows():
                    return False
                if self._cols != other.cols():
                    return False
                if isinstance(other, BitMat):
                    m2 = other._mat
                    for i in range(self._rows):
                        for j in range(self._bcols):
                            if self._mat[i][j] != m2[i][j]:
                                return False
                    return True
                else:
                    for i in range(self._rows):
                        if self.get_row(i) != other.get_row(i):
                            return False
                    return True
        return ConstBaseMatrix.__eq__(self,other)

    def __neq__(self,other):
        mat_class = other.mat_class()
        if isinstance(mat_class, BitMatClass):
            if self._bs == mat_class._bs:
                if self._rows != other.rows():
                    return True
                if self._cols != other.cols():
                    return True
                if isinstance(other, BitMat):
                    m2 = other._mat
                    for i in range(self._rows):
                        for j in range(self._bcols):
                            if self._mat[i][j] != m2[i][j]:
                                return True
                    return False
                else:
                    for i in range(self._rows):
                        if self.get_row(i) != other.get_row(i):
                            return True
                    return False
        return ConstBaseMatrix.__neq__(self,other)

    def is_zero(self):
        for row in self._mat:
            for r in row:
                if r != 0:
                    return False
        return True

    def __repr__(self):
        rows = self._rows
        cols = self._cols
        res = "%d-by-%d bit matrix:\n" % (rows, cols)
        for i in range(rows):
            res += "".join(str(self.get_entry(i,j)) for j in range(cols))
            if i < rows:
                res += "\n"
        return res

    # extra functions

    def rand(rows, cols, block_size=60):
        blocks = cols // block_size
        of = cols % block_size
        res = [[randint(0, (1<<block_size) - 1) for j in range(blocks)] for i in range(rows)]
        if of>0:
            for r in res:
                r.append(randint(0, (1<<of) - 1))
        return BitMat(rows,cols, block_size = block_size, bit_mat = res)


class BitMatClass(BaseMatClass):
    """
    A base class of the metaclass for the type of the given matrix.
    One access the matrix type by calling the method mat_class().
    It should be a singleton,
    but I did not have time to implement a proper singleton design pattern, though.
    """

    def __init__(self, block_size):
        self._bs = block_size
    
    def zero_mat(self, rows, cols):
        """
        Returns a rows-by-cols zero matrix of the given matrix class.
        """
        return BitMat(rows, cols, block_size = self._bs)

    def rand(self, rows, cols):
        """
        Returns a random rows-by-cols matrix of the given matrix class.
        """
        return BitMat.rand(rows, cols, self._bs)

    def vec_class(self):
        """
        Returns the vector class with with the matrix class works (or works best).
        """
        return BitVecClass(self._bs)

    def __eq__(self, other):
        return self._bs == other._bs

    # functions you get for free, but you may want reimplement

    def __neq__(self, other):
        return self._bs != other._bs

    def const_zero_mat(self, rows, cols):
        """
        Returns an immutable rows-by-cols zero matrix of the given matrix class.
        """
        return ZeroBitMat(rows, cols, self._bs) 

    def const_id_mat(self,rows):
        return IdBitMat(rows, self._bs)



class ZeroBitMat(ConstBaseMatrix):
    """
    Immutable zero bit matrix.
    """

    def __init__(self, rows, cols, block_size = 60):
        self._rows = rows 
        self._cols = cols 
        self._bs = block_size
    
    # implementation of functions from ConstBaseMatrix

    def rows(self):
        return self._rows

    def cols(self):
        return self._cols

    def get_entry(self,i,j):
        """
        Returns the (i,j) entry.
        """
        return 0

    def get_row(self,i):
        """
        Returns a vector representing the i-th row.
        """
        return BitVec(self._cols)

    def get_col(self,i):
        """
        Returns a vector representing the i-th column.
        """
        return BitVec(self._rows)

    def Gauss_elim_extended_extras(self,extra):
        """
        @extra is another matrix with the same number of rows as @self.
        Performs Gauss elimination on @self.
        The row operations done on @self are also performed on a copy @extra.
        Returns the Gauss elimination of @self, a list of pairs with the leading entries, and the resulting @extra.
        """
        return BitMat(0, self._cols, block_size = self._bs), [], extra.clone()

    def clone(self):
        return ZeroBitMat(self._rows, self._cols, self._bs)

    def mat_class(self):
        """
        Returns the matrix class of the given object.
        Use the matrix class to access various functions which are not instant specific, e.g.,
        constructors of that class.
        (Note: this is not the same as type().)
        """
        return BitMatClass(self._bs)

    def to_mutable(self):
        """
        Returns a self if mutable. Otherwise returns a multable copy.
        """
        return BitMat(self._rows, self._cols, block_size = self._bs)

    # functions for the bit matrices family of classes

    def _bit_array_clone(self):
        """
        Returns a copy of the internal 2-dim array. Not meant for users.
        """
        bcols = (self._cols + self._bs - 1) // self._bs
        return [[0 for c in range(bcols)] for r in range(self._rows)]

    def get_block_size(self):
        """
        Returns the size of each bit block.
        """
        return self._bs

    def _bit_array(self):
        """
        If mutable, returns actual inner bit array.
        If immutable, returns a copy of inner bit array.
        """
        bcols = (self._cols + self._bs - 1) // self._bs
        return [[0 for c in range(bcols)] for r in range(self._rows)]

    # overriding functions from the base class ConstBaseMatrix

    def transpose(self):
        """
        Returns a matrix representing the transpose of @self.
        Output is immutable.
        """
        return ZeroBitMat(self._cols, self._rows, self._bs)

    def __add__(self,other):
        return other.clone()

    def __sub__(self,other):
        return other.clone() # we're in characteristic 2.

    def __neg__(self):
        return ZeroBitMat(self._rows, self._cols, self._bs)

    def __mul__(self,other):
        if self.cols() != other.rows():
            raise "Dimensions are incompatible!"
        if other.is_mutable():
            return BitMat(self._rows, other.cols(), block_size = self._bs)
        else:
            return ZeroBitMat(self._rows, other.cols(), self._bs)

    def __eq__(self,other):   
        if self._rows != other.rows():
            return False
        if self._cols != other.cols():
            return False
        return other.is_zero()

    def __neq__(self,other):
        if self._rows != other.rows():
            return True
        if self._cols != other.cols():
            return True
        return not other.is_zero()

    def is_zero(self):
        return True

    def apply(self,vec):
        """
        Returns @self * @vec where @vec is a vector regarded as a column vector.
        """
        return vec.vec_class().zero_vec(self._rows) # zero vector

    def __call__(self,vec):
        """
        Returns @self * @vec where @vec is a vector regarded as a column vector.
        """
        return vec.vec_class().zero_vec(self._rows) # zero vector

    def __repr__(self):
        rows = self._rows
        cols = self._cols
        line = "0" * rows + "\n"
        res = "%d-by-%d bit matrix:\n" % (rows, cols) + (line * cols)
        return res
    

class IdBitMat(ConstBaseMatrix):
    """
    An immutable identity bit matrix.
    """

    def __init__(self, rows, block_size = 60):
        self._rows = rows 
        self._bs = block_size
    
    # implementation of functions of ConstBaseMatrix

    def rows(self):
        return self._rows

    def cols(self):
        return self._rows

    def get_entry(self,i,j):
        """
        Returns the (i,j) entry.
        """
        if i == j:
            return 1
        return 0

    def get_row(self,i):
        """
        Returns a vector representing the i-th row.
        """
        res = BitVec(self._rows)
        res.set_entry(i,1)
        return res

    def get_col(self,i):
        """
        Returns a vector representing the i-th column.
        """
        res = BitVec(self._rows)
        res.set_entry(i,1)
        return res

    def Gauss_elim_extended_extras(self,extra):
        """
        @extra is another matrix with the same number of rows as @self.
        Performs Gauss elimination on @self.
        The row operations done on @self are also performed on a copy @extra.
        Returns the Gauss elimination of @self, a list of pairs with the leading entries, and the resulting @extra.
        """
        return IdBitMat(self._rows,self._bs), [(i,i) for i in range(self._rows)], extra.clone()

    def clone(self):
        return IdBitMat(self._rows, self._bs)

    def mat_class(self):
        """
        Returns the matrix class of the given object.
        Use the matrix class to access various functions which are not instant specific, e.g.,
        constructors of that class.
        (Note: this is not the same as type().)
        """
        return BitMatClass(self._bs)

    def to_mutable(self):
        """
        Returns a self if mutable. Otherwise returns a multable copy.
        """
        rows = self._rows
        res = BitMat(rows, rows, block_size = self._bs)
        for i in range(rows):
            res[i,i] = 1
        return res

    # functions for the bit matrices family of classes

    def _bit_array_clone(self):
        """
        Returns a copy of the internal 2-dim array. Not meant for users.
        """
        rows = self._rows
        bs = self._bs
        bcols = (rows + bs - 1) // bs
        res = [[0 for c in range(bcols)] for r in range(rows)]
        for i in range(rows):
            res[i][i // bs] |= (1 << (i%bs) )
        return res

    def get_block_size(self):
        """
        Returns the size of each bit block.
        """
        return self._bs

    def _bit_array(self):
        """
        If mutable, returns actual inner bit array.
        If immutable, returns a copy of inner bit array.
        """
        rows = self._rows
        bs = self._bs
        bcols = (rows + bs - 1) // bs
        res = [[0 for c in range(bcols)] for r in range(rows)]
        for i in range(rows):
            res[i][i // bs] |= (1 << (i%bs) )
        return res

    # overriding functions from the base class ConstBaseMatrix

    def transpose(self):
        """
        Returns a matrix representing the transpose of @self.
        Output is immutable.
        """
        return IdBitMat(self._rows, self._bs)

    def __add__(self,other):
        res = other.mutable_clone()
        for i in range(self._rows):
            res[i,i] = 1 - other[i,i]
        return res

    def __sub__(self,other):
        res = other.mutable_clone()
        for i in range(self._rows):
            res[i,i] = 1 - other[i,i]
        return res

    def __neg__(self):
        return IdBitMat(self._rows, self._bs) # we're in characteristic 2

    def __mul__(self,other):
        if self._rows != other.rows():
            raise "Dimensions are incompatible!"
        return other.clone()

    def __eq__(self,other):
        if isinstance(other, BitMat):
            other_rows = other._rows
            if self._rows != other_rows:
                return False
            if self._rows != other._cols:
                return False
            m2 = other._mat
            bs = other._bs
            bcols = (other_rows + bs - 1) // bs
            for i in range(self._rows):
                special = i//bs
                for j in range(bcols):
                    if j != special:
                        if m2[i][j] != 0:
                            return False
                    else:
                        if m2[i][j] != (1 << (i%bs)):
                            return False
            return True
        if isinstance(other, IdBitMat):
            return self._rows == other._rows # we already checked that the number of rows match
        if isinstance(other, ZeroBitMat):
            return self._rows == 0 and other._rows == 0 and other._cols == 0
        return ConstBaseMatrix.__eq__(self, other)

    def __neq__(self,other):
        if isinstance(other, BitMat):
            other_rows = other._rows
            if self._rows != other_rows:
                return True
            if self._rows != other._cols:
                return True
            m2 = other._mat
            bs = other._bs
            bcols = (other_rows + bs - 1) // bs
            for i in range(self._rows):
                special = i//bs
                for j in range(bcols):
                    if j != special:
                        if m2[i][j] != 0:
                            return True
                    else:
                        if m2[i][j] != (1 << (i%bs)):
                            return True
            return False
        if isinstance(other, IdBitMat):
            return self._rows != other._rows # we already checked that the number of rows match
        if isinstance(other, ZeroBitMat):
            return self._rows != 0 or other._rows != 0 or other._cols != 0
        return ConstBaseMatrix.__neq__(self, other)

    def is_zero(self):
        return self._rows == 0

    def apply(self,vec):
        """
        Returns @self * @vec where @vec is a vector regarded as a column vector.
        """
        return vec.clone()

    def __call__(self,vec):
        """
        Returns @self * @vec where @vec is a vector regarded as a column vector.
        """
        return vec.clone()

    def __repr__(self):
        rows = self._rows
        res = "%d-by-%d bit matrix:\n" % (rows, rows)
        for i in range(rows):
            res += "0"*i + "1" + "0"*(rows-i-1) + "\n"
        return res



