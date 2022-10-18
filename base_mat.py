


class BaseVector(object):
    """
    Base class for vector. Do not use unless you really know what you are doing.
    """

    # functions you need to implement

    def get_entry(self,i):
        """
        Returns the i-th entry of the vector.
        """
        pass

    def set_entry(self,i,v):
        """
        Returns the i-th entry to v
        """
        pass

    def __len__(self):
        pass

    def vec_class(self):
        """
        Returns the vector class of the given object.
        Use the vector class to access various functions which are not instant specific, e.g.,
        constructors of that class.
        (Note: this is not the same as type().)
        """
        pass

    def clone(self):
        pass

    def scale(self, scalar):
        """
        Returns the vector obtained from scaling @self by the given scalar.
        """
        pass

    # functions you get for free (although you might want to make them more efficient)

    def __getitem__(self,key):
        return self.get_entry(key)

    def __setitem__(self,key,value):
        self.set_entry(key,value)

    def get_segment(self,i,j):
        """
        Returns the i,...,j-1 indices of the vector as a new vector.
        """
        res = self.vec_class().zero_vec(j-i)
        for k in range(i,j):
            res.set_entry(k-i, self.get_entry(k))
        return res

    def set_segment(self,i,vec):
        """
        Sets the entries from i onward to vec
        """
        for k in range(i,i+len(vec)):
            self.set_entry(k,vec.get_entry(k-i))

    def __add__(self, other):
        l = len(self)
        res = self.vec_class().zero_vec(l)
        for i in range(l):
            res.set_entry(i, self.get_entry(i) + other.get_entry(i))
        return res

    def __sub__(self, other):
        l = len(self)
        res = self.vec_class().zero_vec(l)
        for i in range(l):
            res.set_entry(i, self.get_entry(i) - other.get_entry(i))
        return res

    def __neg__(self):
        l = len(self)
        res = self.vec_class().zero_vec(l)
        for i in range(l):
            res.set_entry(i, -self.get_entry(i))
        return res

    def __eq__(self,other):
        if len(self) != len(other):
            return False
        for i in range(len(self)):
            if self.get_entry(i) != other.get_entry(i):
                return False
        return True

    def __neq__(self,other):
        if len(self) != len(other):
            return True
        for i in range(len(self)):
            if self.get_entry(i) != other.get_entry(i):
                return True
        return False

    def is_zero(self):
        return self == self.vec_class().zero_vec(len(self))

    def __mul__(self,other):
        """
        Returns the dot-product of @self and @other.
        """
        res = self.vec_class().zero_vec(1).get_entry(0) # fishy way of getting zero.
        for i in range(len(self)):
            res = res + self.get_entry(i) * other.get_entry(i)
        return res

    def __repr__(self):
        return "[" + (",".join(str(self.get_entry(i)) for i in range(len(self)))) + "]"

    def __iter__(self):
        """
        Returns an iterable for the entries of the vector.
        """
        for i in range(len(self)):
            yield self.get_entry(i)

class BaseVecClass(object):
    """
    A base class of the metaclass for the type of the BaseVector.
    One access the vector type by calling the method vec_class().
    This class holds methods which are not instant specific.
    Iheriting classes should be a singletons.
    I did not have time to implement a proper singleton design pattern, though.
    """

    def zero_vec(self,n):
        """
        Returns a zero vector of length n.
        """
        pass

    def rand(self,n):
        """
        Returns a random vector of length n.
        """
        pass

    def __eq__(self,other):
        pass

    # functions you get for free:

    def __neq__(self,other):
        return not self==other


def linear_combination(coeffs, vecs, zero_vec=None):
    """
    Returns a linear combination of the given vectors relative to the given coefficients.
    Provide the @zero_vec in case the list might be empty.
    """
    if zero_vec is None:
        res = vecs[0].vec_class().zero_vec(len(vecs[0]))
    else:
        res = zero_vec.clone()
    for i in range(len(coeffs)):
        res = res + vecs[i].scale(coeffs[i])
    return res

def rand_linear_combination(vecs, zero_vec=None):
    """
    Returns a random linear combination of the given vectors.
    Provide the @zero_vec in case the list might be empty.
    """
    if zero_vec is None:
        vec_class = vecs[0].vec_class()
        res = vec_class.zero_vec(len(vecs[0]))
    else:
        vec_class = zero_vec.vec_class()
        res = zero_vec.clone()
    coeffs = vec_class.rand(len(vecs))
    for i in range(len(vecs)):
        res = res + vecs[i].scale(coeffs[i])
    return res


class ConstBaseMatrix(object):
    """
    Base class for an immutable matrix. Do not use unless you really know what you are doing.
    """

    # functions to implement

    def rows(self):
        pass

    def cols(self):
        pass

    def get_entry(self,i,j):
        """
        Returns the (i,j) entry.
        """
        pass

    def Gauss_elim_extended_extras(self,extra):
        """
        @extra is another matrix with the same number of rows as @self.
        Performs Gauss elimination on @self.
        The row operations done on @self are also performed on a copy @extra.
        Returns the Gauss elimination of @self, a list of pairs with the leading entries, and the resulting @extra.
        Output matrices may be immutable if @self is.
        """
        pass

    def clone(self):
        pass

    def mat_class(self):
        """
        Returns the matrix class of the given object.
        Use the matrix class to access various functions which are not instant specific, e.g.,
        constructors of that class.
        (Note: this is not the same as type().)
        """
        pass

    def to_mutable(self):
        """
        Returns a self if mutable. Otherwise returns a multable copy.
        """
        pass

    # functions you get for free if you implemented all of the above (although you may want to make them more efficient).

    def __getitem__(self,key):
        i,j = key
        return self.get_entry(i,j)

    def get_row(self,i):
        """
        Returns a vector representing the i-th row.
        """
        cols = self.cols()
        res = self.mat_class().vec_class().zero_vec(cols)
        for j in range(cols):
            res.set_entry(j, self.get_entry(i,j))
        return res

    def get_rows(self):
        """
        Returns an iterable enumerating vectors representing the rows of the matrix.
        """
        for i in range(self.rows()):
            yield self.get_row(i)

    def get_row_as(self, i, vec_class):
        """
        Returns a vector representing the i-th row. The result has the given vector class.
        """
        cols = self.cols()
        res = vec_class.zero_vec(cols)
        for j in range(cols):
            res.set_entry(j, self.get_entry(i,j))
        return res

    def get_col(self,i):
        """
        Returns a vector representing the i-th column.
        """
        rows = self.rows()
        res = self.mat_class().vec_class().zero_vec(rows)
        for j in range(rows):
            res.set_entry(j, self.get_entry(j,i))
        return res

    def get_cols(self):
        """
        Returns an iterable enumerating vectors representing the columns of the matrix.
        """
        for i in range(self.cols()):
            yield self.get_col(i)

    def get_col_as(self, i, vec_class):
        """
        Returns a vector representing the i-th column. The result has the given vector class.
        """
        rows = self.rows()
        res = vec_class.zero_vec(rows)
        for j in range(rows):
            res.set_entry(j, self.get_entry(j,i))
        return res

    def get_rect(self, i1, i2, j1, j2):
        """
        Returns a matrix representing the submatrix consisting of the (i,j)-entries
        with i1 <= i < i2 and j1 <= j < j2. 
        """
        res = self.mat_class().zero_mat(i2-i1, j2-j1)
        for i in range(i1, i2):
            for j in range(j1, j2):
                res.set_entry(i-i1, j-j1, self.get_entry(i,j))
        return res

    def transpose(self):
        """
        Returns a matrix representing the transpose of @self.
        Output matrix may be immutable if @self is.
        """
        rows = self.rows()
        cols = self.cols()
        res = self.mat_class().zero_mat(cols, rows)
        for i in range(rows):
            for j in range(cols):
                res.set_entry(j,i,self.get_entry(i,j))
        return res

    def Gauss_elim_extended(self):
        """
        Returns the Gauss elimination of the matrix and a list of the leading entries (pairs of coordinates)
        Output matrix may be immutable if @self is.
        """
        extra = self.mat_class().const_zero_mat(self.rows(),0)
        res, leading, E = self.Gauss_elim_extended_extras(extra)
        return res, leading

    def Gauss_elim(self):
        """
        Returns the Gauss elimination of the matrix.
        Output matrix may be immutable if @self is.
        """
        res, leading = self.Gauss_elim_extended()
        return res

    def null_space(self):
        """
        Returns a matrix representing a basis to the (right) null space of a matrix M.
        Output matrix may be immutable if @self is.
        """
        E, leading = self.Gauss_elim_extended()
        cols = self.cols()
        res = self.mat_class().zero_mat(cols - len(leading), cols)
        L = {}
        for r,c in leading:
            L[c] = r
        j=0
        for i in range(cols):
            if i not in L:
                res.set_to_one(j,i)
                for k,l in L.items():
                    res.set_entry(j,k, -E.get_entry(l,i))
                j += 1
        return res

    def quasi_inverse(self):
        """
        Returns a quasi-inverse of @self, i.e., a matrix N such that self*N*self = self.
        If a system M*x = b is known to have a solution, then x = N*b is one such solution.
        The qausi-inverse is also a left / right inverse of @self, if such an inverse exists.
        Output matrix may be immutable if @self is.
        """
        rows = self.rows()
        I = self.mat_class().const_id_mat(rows)
        G,leading,P = self.Gauss_elim_extended_extras(I)
        Q = self.mat_class().zero_mat(self.cols(), rows) 
        for i,j in leading:
            Q.set_row(j,P.get_row(i))
            #for k in range(rows)
            #    Q[j,k] = P.[i,k]
        return Q

    def inverse(self):
        """
        Returns the inverse of the matrix or None if the matrix is not invertible.
        """
        rows = self.rows()
        if rows != self.cols():
            return None
        I = self.mat_class().const_id_mat(rows)
        G,leading,P = self.Gauss_elim_extended_extras(I)
        if G.rows() != rows:
            return None
        return P
        
    def ker(self):
        """
        Returns a matrix which represents a kernel of @self.
        That is @self * res = 0, res is one-to-one as a linear transformation, and image(res) = null space of @self.
        Output matrix may be immutable if @self is.
        """
        return self.null_space().transpose()

    def coker(self):
        """
        Returns a matrix which represents a cokernel of @self.
        That is res * @self = 0, res is onto as a linear transformation, and ker(res) = column space of @self.
        Output matrix may be immutable if @self is.
        """
        return self.transpose().null_space()

    def __add__(self,other):
        rows = self.rows()
        cols = self.cols()
        res = self.mat_class().zero_mat(rows, cols)
        for i in range(rows):
            for j in range(cols):
                res.set_entry(i,j, self.get_entry(i,j) + other.get_entry(i,j))
        return res

    def __sub__(self,other):
        rows = self.rows()
        cols = self.cols()
        res = self.mat_class().zero_mat(rows, cols)
        for i in range(rows):
            for j in range(cols):
                res.set_entry(i,j, self.get_entry(i,j) - other.get_entry(i,j))
        return res

    def __neg__(self):
        rows = self.rows()
        cols = self.cols()
        res = self.mat_class().zero_mat(rows, cols)
        for i in range(rows):
            for j in range(cols):
                res.set_entry(i,j, -self.get_entry(i,j))
        return res

    def __mul__(self,other):
        if self.cols() != other.rows():
            raise "Dimensions are incompatible!"
        rows = self.rows()
        cols = other.cols()
        vec_class = self.mat_class().vec_class()
        res = self.mat_class().zero_mat(rows, cols)
        for j in range(cols):
            col = other.get_col_as(j, vec_class)
            for i in range(rows):
                res.set_entry(i,j, self.get_row(i) * col)
        return res

    def __eq__(self,other):
        if self.rows() != other.rows():
            return False
        if self.cols() != other.cols():
            return False
        for i in range(self.rows()):
            for j in range(self.cols()):
                if self.get_entry(i,j) != other.get_entry(i,j):
                    return False
        return True

    def __neq__(self,other):
        if self.rows() != other.rows():
            return True
        if self.cols() != other.cols():
            return True
        for i in range(self.rows()):
            for j in range(self.cols()):
                if self.get_entry(i,j) != other.get_entry(i,j):
                    return True
        return False

    def is_zero(self):
        return self == self.mat_class().const_zero_mat(self.rows(), self.cols()) 

    def apply(self,vec):
        """
        Returns @self * @vec where @vec is a vector regarded as a column vector.
        """
        rows = self.rows()
        res = vec.vec_class().zero_vec(rows)
        for i in range(rows):
            res.set_entry(i, self.get_row(i)*vec)
        return res

    def __call__(self,vec):
        """
        Returns @self * @vec where @vec is a vector regarded as a column vector.
        """
        return self.apply(vec)

    def __repr__(self):
        rows = self.rows()
        cols = self.cols()
        res = "%d-by-%d matrix:\n" % (rows, cols)
        for i in range(rows):
            res += "[" + (",".join(str(self.get_entry(i,j)) for j in range(cols))) + "]"
            if i < rows:
                res += "\n"
        return res

    def rank(self):
        return self.Gauss_elim().rows()
    
    def is_mutable(self):
        return False

    def mutable_clone(self):
        """
        Returns a mutable copy.
        """
        return self.clone().to_mutable() # this implementation is also good for class inheriting from BaseMatrix

class BaseMatrix(ConstBaseMatrix):
    """
    Base class for a mutable matrix. Do not use unless you really know what you are doing.
    """

    # functions to implement include those in ConstBaseMatrix and the following:

    def set_entry(self,i,j,v):
        """
        Sets entry (i,j) to v.
        """

    def set_to_zero(self,i,j):
        """
        Sets entry (i,j) to 0.
        """
        pass

    def set_to_one(self,i,j):
        """
        Sets entry (i,j) to 1.
        """
        pass

    def add_zero_rows(self, rows):
        """
        Adds @rows rows of zeroes to @self.
        """
        pass

    def append_rows(self, other):
        """
        Adds the rows of @other at the end of @self.
        """
        pass

    # functions you get for free if you implemented all of the above (although you may want to make them more efficient).

    def to_mutable(self):
        """
        Returns a self if mutable. Otherwise returns a multable copy.
        """
        return self

    def __setitem__(self,key,value):
        i,j = key
        self.set_entry(i,j,value)

    def set_rect(self,i,j,mat):
        """
        Sets the rectangle strating at (i,j) and having @mat.rows() rows and @mat.cols() columns to @mat.
        """
        for r in range(mat.rows()):
            for c in range(mat.cols()):
                self.set_entry(i+r, j+c, mat.get_entry(r,c))

    def set_row(self,i,vec):
        """
        Sets the i-th row to vec.
        """
        for j in range(len(vec)):
            self.set_entry(i,j, vec.get_entry(j))

    def set_col(self,i,vec):
        """
        Sets the i-th column to vec.
        """
        for j in range(len(vec)):
            self.set_entry(j,i, vec.get_entry(j))

    def set_row_segment(self,i,j,vec):
        """
        Sets the entries starting at (i,j) and to its right to @vec.
        """
        for c in range(len(vec)):
            self.set_entry(i, j+c, vec.get_entry(c))

    def set_col_segment(self,i,j,vec):
        """
        Sets the entries starting at (i,j) and below to @vec.
        """
        for c in range(len(vec)):
            self.set_entry(i+c, j, vec.get_entry(c))

    def is_mutable(self):
        return True

class BaseMatClass(object):
    """
    A base class of the metaclass for the type of the given matrix.
    One access the matrix type by calling the method mat_class().
    Iheriting classes should be a singletons.
    I did not have time to implement a proper singleton design pattern, though.
    """
    
    def zero_mat(self, rows, cols):
        """
        Returns a rows-by-cols zero matrix of the given matrix class.
        """
        pass

    def rand(self, rows, cols):
        """
        Returns a random rows-by-cols matrix of the given matrix class.
        """
        pass

    def vec_class(self):
        """
        Returns the vector class with with the matrix class works (or works best).
        """
        pass

    def __eq__(self, other):
        pass

    # functions you get for free, but you may want reimplement

    def __neq__(self, other):
        return not self==other

    def const_zero_mat(self, rows, cols):
        """
        Returns an immutable rows-by-cols zero matrix of the given matrix class.
        """
        return self.zero_mat(rows, cols)

    def id_mat(self,rows):
        """
        Returns a mutable identity matrix of the given dimensions and of the same type as @self.
        """
        res = self.zero_mat(rows,rows)
        for i in range(rows):
            res.set_to_one(i,i)
        return res

    def const_id_mat(self,rows):
        return self.id_mat(rows)

    def rand_inv_mat(self, rows):
        res = self.rand(rows,rows)
        while res.rank() != rows:
            res = self.rand(rows,rows)
        return res
    
    def cohomology_basis(self,d1,d0):
        """
        Receives two matrices d0, d1 such that d1 * d0 = 0.
        Returns a matrix whose rows span a complement of im(d0) inside ker(d1).
        """
        B,leading = d0.transpose().Gauss_elim_extended()
        del B
        D1 = d1.mutable_clone()
        D1.add_zero_rows(len(leading))
        r = d1.rows()
        for k,i in leading:
            for j in range(r):
                D1.set_to_zero(j,i)
            D1.set_to_one(r + k, i)
        return D1.null_space()

    def cohomology_dual_basis(self,d1,d0):
        """
        Receives two matrices d0, d1 such that d1 * d0 = 0.
        Returns a matrix such that its rows v1,..,vt define an isomorphism
        x -> (<v1,x>,...,<vt,x>) from ker(d1)/im(d0) to (field)^t.
        """
        B, leading = d1.Gauss_elim_extended()
        del B
        d0t = d0.transpose().to_mutable()
        d0t.add_zero_rows(len(leading))
        c = d0.cols()
        for k,i in leading:
            for j in range(c):
                d0t.set_to_zero(j,i)
            d0t.set_to_one(c + k, i)
        return d0t.null_space()

    def rand_cohomology_basis(self,d1,d0):
        """
        Receives two matrices d0, d1 such that d1 * d0 = 0.
        Returns a random matrix whose rows span a complement of im(d0) inside ker(d1).
        """
        B,leading = d0.transpose().Gauss_elim_extended()
        D1 = d1.mutable_clone()
        D1.add_zero_rows(len(leading))
        r = d1.rows()
        for k,i in leading:
            for j in range(r):
                D1.set_to_zero(j,i)
            D1.set_to_one(r + k, i)
        C = D1.null_space()
        rows_C = C.rows()
        C.append_rows(B)
        del B
        del D1
        E = self.rand_inv_mat(rows_C)
        M = self.rand(rows_C, C.rows())
        for i in range(E.rows()):
            M.set_row_segment(i,0, E.get_row(i))
        return M * C


def cohomology_basis(d1,d0):
    """
    Receives two matrices d0, d1 such that d1 * d0 = 0.
    Returns a matrix whose rows span a complement of im(d0) inside ker(d1).
    """
    return d1.mat_class().cohomology_basis(d1,d0)

def cohomology_dual_basis(d1,d0):
    """
    Receives two matrices d0, d1 such that d1 * d0 = 0.
    Returns a matrix such that its rows v1,..,vt define an isomorphism
    x -> (<v1,x>,...,<vt,x>) from ker(d1)/im(d0) to (field)^t.
    """
    return d1.mat_class().cohomology_dual_basis(d1,d0)

def rand_cohomology_basis(d1,d0):
    """
    Receives two matrices d0, d1 such that d1 * d0 = 0.
    Returns a random matrix whose rows span a complement of im(d0) inside ker(d1).
    """
    return d1.mat_class().rand_cohomology_basis(d1,d0)



    
