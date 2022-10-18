from base_mat import *


class Blocks(object):
    """
    Suppose we partition [0,...,N-1] into interval and name each interval with a key.
    This class represents this data structure, allowing for conversion between keys and subintervals.
    """

    def __init__(self, block_lens, key_order=None):
        """
        @block_lens is an iterable of pairs (key,lenght) consisting of a key and the length of the block it represents.
        @key_order is the an optional iterable specifying an order of the keys. If left unspecified,
            an arbitrary order is taken.
        """
        if key_order is None:
            self._keys = []
        else:
            self._keys = [k for k in key_order]
        self._block_lens = {}
        for k,v in block_lens:
            self._block_lens[k] = v
            if key_order is None:
                self._keys.append(k)
        self._intervals = {}
        self._index_to_key = []
        of = 0
        for k in self._keys:
            next_of = of + self._block_lens[k]
            self._intervals[k] = (of,next_of)
            for i in range(0,next_of-of):
                self._index_to_key.append((k,i))
            of = next_of
        self._len = of

    def __eq__(self,other):
        if self is other:
            return True # saving some time
        return self._keys == other._keys and self._block_lens == other._block_lens

    def __neq__(self,other):
        if self is other:
            return False # saving some time
        return self._keys != other._keys or self._block_lens != other._block_lens

    def key_to_index(self,key,j=0):
        """
        Returns the index of the j-th posion in the key-th block.
        """
        return self._intervals[key][0]+j

    def index_to_key(self,i):
        """
        Returns key,j for which the j-th position in the key-th block corresponds to index i.
        """
        return self._index_to_key[i]

    def interval(self,key):
        return self._intervals[key]

    def __len__(self):
        """
        The total length of all the blocks.
        """
        return self._len

    def block_len(self,key):
        """
        The length of of the key-th block
        """
        return self._block_lens[key]

    def keys(self):
        """
        Returns an iterable giving all keys in order.
        """
        for key in self._keys:
            yield key

class BlockVec(BaseVector):
    """
    Represents a block vector. The keys can be any hashable.
    Note:
    1. A BlockVector is a an ordinary vector, and thus supports all operations on vectors
        supplied by the BaseVector class. However, the operations such as addition will return
        non-blocked vectors if the one of the arguments admits no block structure, or
        if the block structures do not match.
    2. The vector class of a block vector is the vector class of its unlying vector.
    """
    
    def __init__(self, blocks, vec=None, vec_class=None):
        """
        @blocks is an instance of Blocks
        @vec is an instance of BaseVector; if it is a BlockVec, then the block structure is stripped off.
        @vec_block the vector class of vec. This variable is mandatory if @vec is None.
        If @vec is omitted and vec_class is given, then a zero vector of the given class is created.
        """
        while isinstance(vec,BlockVec):
            vec = vec._vec
        self._vec = vec
        if vec is None:
            self._vec = vec_class.zero_vec(len(blocks))
        self._blocks = blocks

    # functions you need to implement

    def get_entry(self,i):
        """
        Returns the i-th entry of the vector.
        """
        return self._vec.get_entry(i)

    def set_entry(self,i,v):
        """
        Returns the i-th entry to v
        """
        self._vec.set_entry(i,v)

    def __len__(self):
        return len(self._vec)

    def vec_class(self):
        """
        Returns the vector class of the given object.
        Use the vector class to access various functions which are not instant specific, e.g.,
        constructors of that class.
        (Note: this is not the same as type().)
        """
        return self._vec.vec_class()

    def clone(self):
        return BlockVec(self._blocks, self._vec)

    def scale(self, scalar):
        """
        Returns the vector obtained from scaling @self by the given scalar.
        """
        return BlockVec(self._blocks, self._vec.scale(scalar))

    #overriding some BaseVector functions

    def __getitem__(self,key):
        return self._vec.get_entry(key)

    def __setitem__(self,key,value):
        self._vec.set_entry(key,value)

    def get_segment(self,i,j):
        """
        Returns the i,...,j-1 indices of the vector as a new vector.
        """
        return self._vec.get_segment(i,j)

    def set_segment(self,i,vec):
        """
        Sets the entries from i onward to vec
        """
        self._vec.set_segment(i,vec)

    def __add__(self, other):
        if not isinstance(other, BlockVec):
            return self._vec + other # other is not a block vector, so block structure is lost
        if self._blocks == other._blocks:
            return BlockVec(self._blocks, self._vec + other._vec)
        return self._vec + other._vec # block structures do not match, so block structure is lost

    def __sub__(self, other):
        if not isinstance(other, BlockVec):
            return self._vec - other # other is not a block vector, so block structure is lost
        if self._blocks == other._blocks:
            return BlockVec(self._blocks, self._vec - other._vec)
        return self._vec - other._vec # block structures do not match, so block structure is lost

    def __neg__(self):
        return BlockVec(self._blocks, -self._vec)

    def __eq__(self,other):
        """
        Equates the two vectors ignoring the block structure.
        """
        if isinstance(other, BlockVec):
            return self._vec == other._vec
        return self._vec == other 

    def __neq__(self,other):
        """
        Equates the two vectors ignoring the block structure.
        """
        if isinstance(other, BlockVec):
            return self._vec != other._vec
        return self._vec != other 

    def is_zero(self):
        return self._vec.is_zero()

    def __mul__(self,other):
        """
        Returns the dot-product of @self and @other.
        """
        if isinstance(other, BlockVec):
            return self._vec * other._vec
        return self._vec * other 

    def __repr__(self):
        b = self._blocks
        res = "["
        for key in b.keys():
            i,j = b.interval(key)
            res += str(key) + ":" + self._vec.get_segment(i,j).__repr__() + ", "
        if len(b) > 0:
            res = res[:-2]
        res += "]"
        return res
    
    # extra functions

    def get_block(self, key):
        i,j = self._blocks.interval(key)
        return self._vec.get_segment(i,j)

    def set_block(self, key, vec):
        i,j = self._blocks.interval(key)
        self._vec.set_segment(i,vec)

    def __call__(self, key):
        """
        Returns a vector reprsenting to block corresponding to the given key.
        Same as get_block.
        """
        i,j = self._blocks.interval(key)
        return self._vec.get_segment(i,j)




class SparseBlockMat(ConstBaseMatrix):
    """
    Represents an immutable sparse block matrix. The keys can be any hashable.
    """
    
    def __init__(self, mat_blocks, row_blocks, col_blocks, mat_class = None):
        """
        @mat_blocks is a dictionary mapping a pair (row_key,col_key) to a BitMat object.
            Missing blocks are interpreted as zero blocks.
        @row_blocks is a Blocks object representing the row blocks
        @col_blocks is a Blocks object representing the column blocks
        @mat_class is the relevant matrix class. If not specified, then it read from some block in mat_blocks
        """
        self._blocks = mat_blocks
        self._row_blocks = row_blocks
        self._col_blocks = col_blocks
        if mat_class is None:
            for b in mat_blocks.values():
                self._mat_class = b.mat_class()
                break
        else:
            self._mat_class = mat_class

    # methods from ConstBaseMatrix we must implement

    def rows(self):
        return len(self._row_blocks)

    def cols(self):
        return len(self._col_blocks)

    def get_entry(self,i,j):
        """
        Return s the (i,j) entry.
        """
        rk, a = self._row_blocks.index_to_key(i)
        ck, b = self._col_blocks.index_to_key(j)
        if (rk, ck) in self._blocks:
            return self._blocks[rk,ck].get_entry(a,b)
        else:
            return self._mat_class.zero_mat(1,1).get_entry(0,0) # a strage way of getting zero

    def get_row(self,i):
        """
        Returns a vector representing the i-th row.
        """
        rk, j = self._row_blocks.index_to_key(i)
        res = self._mat_class.vec_class().zero_vec(len(self._col_blocks))
        for ck in self._col_blocks.keys():
            if (rk,ck) in self._blocks:
                u,v = self._col_blocks.interval(ck)
                res.set_segment(u, self._blocks[rk,ck].get_row(j))
        return res

    def get_col(self,i):
        """
        Returns a vector representing the i-th column.
        """
        ck, j = self._col_blocks.index_to_key(i)
        res = self._mat_class.vec_class().zero_vec(len(self._row_blocks))
        for rk in self._row_blocks.keys():
            if (rk,ck) in self._blocks:
                u,v = self._row_blocks.interval(rk)
                res.set_segment(u, self._blocks[rk,ck].get_col(j))
        return res

    def Gauss_elim_extended_extras(self,extra):
        """
        @extra is another matrix with the same number of rows as @self.
        Performs Gauss elimination on @self.
        The row operations done on @self are also performed on a copy @extra.
        Returns the Gauss elimination of @self, a list of pairs with the leading entries, and the resulting @extra.
        Output matrices may be immutable if @self is.
        """
        if isinstance(extra, SparseBlockMat):
            return self.forget_blocks().Gauss_elim_extended_extras(extra.forget_blocks())
        return self.forget_blocks().Gauss_elim_extended_extras(extra)

    def clone(self):
        return SparseBlockMat(self._blocks, self._row_blocks, self._col_blocks, self._mat_class)

    def mat_class(self):
        """
        Returns the matrix class of the given object.
        Use the matrix class to access various functions which are not instant specific, e.g.,
        constructors of that class.
        (Note: this is not the same as type().)
        """
        return self._mat_class

    def to_mutable(self):
        """
        Returns a self if mutable. Otherwise returns a multable copy.
        """
        return self.forget_blocks()

    # methods overring those of ConstBaseMatrix

    def Gauss_elim_extended(self):
        """
        Returns the Gauss elimination of the matrix and a list of the leading entries (pairs of coordinates)
        Output matrix may be immutable if @self is.
        """
        return self.forget_blocks().Gauss_elim_extended()

    def Gauss_elim(self):
        """
        Returns the Gauss elimination of the matrix.
        Output matrix may be immutable if @self is.
        """
        return self.forget_blocks().Gauss_elim()

    def __add__(self,other):
        if self.rows() != other.rows() or self.cols() != other.cols():
            raise "Dimensions are incompatible!"
        if not isinstance(other, SparseBlockMat):
            return self.forget_blocks() + other # other is not a block matrix, so block structure is lost.
        if (self._row_blocks != other._row_blocks) or (self._col_blocks != other._col_blocks):
            return self.forget_blocks() + other.forget_blocks() # blocks do not match, so block structure is lost.
        # block structures match
        res_blocks = {k:b.clone() for k,b in self._blocks.items()} 
        for k,b in other._blocks.items():
            if k in res_blocks:
                res_blocks[k] = res_blocks[k] + other._blocks[k]
            else:
                res_blocks[k] = other._blocks[k]
        return SparseBlockMat(res_blocks, self._row_blocks, self._col_blocks, self._mat_class)
            

    def __sub__(self,other):
        if self.rows() != other.rows() or self.cols() != other.cols():
            raise "Dimensions are incompatible!"
        if not isinstance(other, SparseBlockMat):
            return self.forget_blocks() - other # other is not a block matrix, so block structure is lost.
        if (self._row_blocks != other._row_blocks) or (self._col_blocks != other._col_blocks):
            return self.forget_blocks() - other.forget_blocks() # blocks do not match, so block structure is lost.
        # block structures match
        res_blocks = {k:b.clone() for k,b in self._blocks.items()} 
        for k,b in other._blocks.items():
            if k in res_blocks:
                res_blocks[k] = res_blocks[k] - other._blocks[k]
            else:
                res_blocks[k] = -other._blocks[k]
        return SparseBlockMat(res_blocks, self._row_blocks, self._col_blocks, self._mat_class)

    def __neg__(self):
        res_blocks = {}
        for k,b in self._blocks.items():
            res_blocks[k] = -b
        return SparseBlockMat(res_blocks, self._row_blocks, self._col_blocks, self._mat_class)

    def __mul__(self,other):
        if self.cols() != other.rows():
            raise "Dimensions are incompatible!"
        if not isinstance(other, SparseBlockMat):
            return self.forget_blocks() * other # other is not a block matrix, so block structure is lost.
        if self._col_blocks != other._row_blocks:
            return self.forget_blocks() * other.forget_blocks() # blocks do not match, so block structure is lost.
        # block structures match
        other_blocks_by_col = {}
        for k,b in other._blocks.items():
            rk, ck = k
            if ck in other_blocks_by_col:
                other_blocks_by_col[ck][rk] = b
            else:
                other_blocks_by_col[ck] = {rk:b}
        res_blocks = {}
        for k,b in self._blocks.items():
            rk, ck = k
            for o_ck, o_ck_blocks in other_blocks_by_col.items():
                if ck in o_ck_blocks:
                    if (rk, o_ck) in res_blocks:
                        res_blocks[rk, o_ck] = res_blocks[rk, o_ck] + b*o_ck_blocks[ck]
                    else:
                        res_blocks[rk, o_ck] = b*o_ck_blocks[ck]
        return SparseBlockMat(res_blocks, self._row_blocks, other._col_blocks, self._mat_class)
    
    def __eq__(self,other):
        if self.rows() != other.rows() or self.cols() != other.cols():
            return False
        if not isinstance(other, SparseBlockMat):
            return BaseMatrix.__eq__(self, other)
        if (self._row_blocks != other._row_blocks) or (self._col_blocks != other._col_blocks):
            return BaseMatrix.__eq__(self, other) # blocks do not match, so block structure is lost.
        # block structures match
        for k,b in self._blocks.items():
            rk, ck = k
            if b != other.get_block(rk,ck):
                return False
        for k,b in other._blocks.items():
            rk, ck = k
            if b != self.get_block(rk,ck):
                return False
        return True
    
    def __neq__(self,other):
        if self.rows() != other.rows() or self.cols() != other.cols():
            return True
        if not isinstance(other, SparseBlockMat):
            return BaseMatrix.__neq__(self, other)
        if (self._row_blocks != other._row_blocks) or (self._col_blocks != other._col_blocks):
            return BaseMatrix.__neq__(self, other) # blocks do not match, so block structure is lost.
        # block structures match
        for k,b in self._blocks.items():
            rk, ck = k
            if b != other.get_block(rk,ck):
                return True
        for k,b in other._blocks:
            rk, ck = k
            if b != self.get_block(rk,ck):
                return True
        return False

    def is_zero(self):
        for b in self._blocks.values():
            if not b.is_zero():
                return False
        return True

    def apply(self,vec):
        """
        Returns @self * @vec where @vec is a vector regarded as a column vector.
        The result is a vector of the same vector class as the matrix. NOT a block vector.
        """
        res = self._mat_class.vec_class().zero_vec(len(self._row_blocks))
        for k,b in self._blocks.items():
            rk,ck = k
            i,j = self._col_blocks.interval(ck)
            k,l = self._row_blocks.interval(rk)
            res.set_segment(k, res.get_segment(k,l) + b.apply(vec.get_segment(i,j)))
        return res

    def __repr__(self):
        return ("%d-by-%d block matrix, blocks =" % (self.rows(), self.cols())) + str(self._blocks)
    
    def is_mutable(self):
        return False


    # additional methods

    def apply_and_return_as_BlockVec(self, vec):
        """
        Returns @self * @vec where @vec is a vector regarded as a column vector.
        The result is a block vector.
        """
        res = BlockVec(blocks = self._row_blocks, vec_class = vec.vec_class())
        for k,b in self._blocks.items():
            rk,ck = k
            i,j = self._col_blocks.interval(ck)
            res.set_block(rk, res.get_block(rk) + b.apply(vec.get_segment(i,j)))
        return res

    def get_block(self,rk,ck):
        """
        Returns the (@rk,@ck)-block.
        """
        if (rk, ck) in self._blocks:
            return self._blocks[rk,ck]
        else:
            return self._mat_class.zero_mat(self._row_blocks.block_len(rk), self._col_blocks.block_len(ck)) # zero block

    def validate(self):
        """
        Makes sure that the blocks are in the right dimensions.
        """
        for k,b in self._blocks.items():
            rk, ck = k
            if b.rows() != self._row_blocks.block_len(rk):
                self._validation_failure = (rk,ck)
                return False
            if b.cols() != self._col_blocks.block_len(ck):
                self._validation_failure = (rk,ck)
                return False
        return True

    def forget_blocks(self):
        """
        Returns a mutable non-block matrix of the relevant matrix class representing @self.
        """
        res = self._mat_class.zero_mat(self.rows(), self.cols())
        for k,b in self._blocks.items():
            rk,ck = k
            r_of = self._row_blocks.key_to_index(rk)
            c_of = self._col_blocks.key_to_index(ck)
            res.set_rect(r_of,c_of,b)
        return res



        

