
# should be used together with the fields module

from base_mat import *
from itertools import chain

class Vector(BaseVector):
    """
    A vector over a generic field.
    """

    def __init__(self, field, size=None, vec=None):
        """
        @field is the class of the field over which the vector is defined.
        @size the the size of the vector
        @vec is a list of the vector coordinates.
        One of @size of @vec should be specified.
        If @vec is omitted, a zero vector is created.
        """
        self._field = field
        if vec is None:
            self._vec = [field.zero()] * size
        else:
            self._vec = [vi for vi in vec]

    # must implement for BaseVector

    def get_entry(self,i):
        """
        Returns the i-th entry of the vector.
        """
        return self._vec[i]

    def set_entry(self,i,v):
        """
        Returns the i-th entry to v
        """
        self._vec[i] = v

    def __len__(self):
        return len(self._vec)

    def vec_class(self):
        """
        Returns the vector class of the given object.
        Use the vector class to access various functions which are not instant specific, e.g.,
        constructors of that class.
        (Note: this is not the same as type().)
        """
        return VectorClass(self._field)

    def clone(self):
        return Vector(self._field, vec = self._vec)

    def scale(self, scalar):
        """
        Returns the vector obtained from scaling @self by the given scalar.
        """
        return Vector(self._field, vec = [scalar*vi for vi in self._vec])

    # additional functions

    def get_field(self):
        return self._field

class VectorClass(BaseVecClass):
    """
    A class representing the vector class of a Vector instance.
    """

    def __init__(self,field):
        self._field = field

    # must implement for BaseVecClass
    
    def zero_vec(self,n):
        """
        Returns a zero vector of length n.
        """
        return Vector(self._field, size = n)

    def rand(self,n):
        """
        Returns a random vector of length n.
        """
        return Vector(self._field, vec = [self._field.rand() for i in range(n)])

    def __eq__(self,other):
        return self._field is other._field

    # overriding functions from BaseVecClass 

    def __neq__(self,other):
        return self._field is not other._field


class Matrix(BaseMatrix):
    """
    A matrix over a generic field.
    """

    def __init__(self, field, rows, cols, mat=None):
        """
        @field is the base field,
        @rows is the number of rows,
        @cols is the number of columns,
        @mat is a list of lists representing the matrix.
        If @mat is omitted, a zero matrix is created.
        """
        self._field = field
        self._cols = cols
        self._rows = rows
        if mat is None:
            zero = field.zero()
            self._mat = [[zero for i in range(cols)] for j in range(rows)]
        else:
            self._mat = [r[:] for r in mat]

    # must implement for ConstBaseMatrix

    def rows(self):
        return self._rows

    def cols(self):
        return self._cols

    def get_entry(self,i,j):
        """
        Returns the (i,j) entry.
        """
        return self._mat[i][j]

    def get_row(self,i):
        """
        Returns a vector representing the i-th row.
        """
        return Vector(self._field, vec=self._mat[i])

    def get_col(self,i):
        """
        Returns a vector representing the i-th column.
        """
        return Vector(self._field, vec=[r[i] for r in self._mat])

    def Gauss_elim_extended_extras(self,extra):
        """
        @extra is another instance of Matrix with the same number of rows as @self.
        Performs Gauss elimination on @self.
        The row operations done on @self are also performed on a copy @extra.
        Returns the Gauss elimination of @self, a list of pairs with the leading entries, and the resulting @extra.
        Output matrices may be immutable if @self is.
        """
        assert self.rows() == extra.rows()
        assert isinstance(extra,Matrix)
        rows = self._rows
        cols = self._cols
        field = self._field
        zero = field.zero()
        mat = self._mat
        if (rows == 0) or (cols == 0):
            return Matrix(self._field, 0, self._cols), [], extra.clone()
        E = [r[:] for r in mat]
        extra_E = [r[:] for r in extra._mat]
        extra_cols = extra.cols()
        leading = []
        i = 0 # leading column
        j = 0 # leading row
        while i < cols:
            s = j
            while s < rows and E[s][i] == zero:
                s += 1
            if s < rows:
                leading.append((j,i))
                E[s], E[j] = E[j], E[s]
                extra_E[s], extra_E[j] = extra_E[j], extra_E[s]
                u = E[j][i].inverse()
                E[j][i] = field.one()
                for k in range(i+1,cols):
                    E[j][k] = E[j][k] * u
                for k in range(extra_cols):
                    extra_E[j][k] = extra_E[j][k] * u
                for l in chain(range(0,j),range(j+1,rows)):
                    if E[l][i] == field(0):
                        continue
                    for k in range(i+1,cols):
                        E[l][k] = E[l][k] - E[l][i] * E[j][k]
                    for k in range(extra_cols):
                        extra_E[l][k] = extra_E[l][k] - E[l][i] * extra_E[j][k]
                    E[l][i] = field(0)
                j = j+1
            i = i+1
        return Matrix(field, j, cols, mat=E[:j]), \
               leading, \
               Matrix(field, len(extra_E), extra_cols, mat=extra_E) # truncate zero rows

    def clone(self):
        return Matrix(self._field, self._rows, self._cols, mat=self._mat)

    def mat_class(self):
        """
        Returns the matrix class of the given object.
        Use the matrix class to access various functions which are not instant specific, e.g.,
        constructors of that class.
        (Note: this is not the same as type().)
        """
        return MatrixClass(self._field)

    def to_mutable(self):
        """
        Returns a self if mutable. Otherwise returns a multable copy.
        """
        return self

    # must implement for BaseMatrix

    def set_entry(self,i,j,v):
        """
        Sets entry (i,j) to v.
        """
        self._mat[i][j] = v

    def set_to_zero(self,i,j):
        """
        Sets entry (i,j) to 0.
        """
        self._mat[i][j] = self._field.zero()

    def set_to_one(self,i,j):
        """
        Sets entry (i,j) to 1.
        """
        self._mat[i][j] = self._field.one()

    def add_zero_rows(self, rows):
        """
        Adds @rows rows of zeroes to @self.
        """
        zero = self._field.zero()
        self._mat += [[zero]*self._cols for i in range(rows)]
        self._rows += rows

    def append_rows(self, other):
        """
        Adds the rows of @other at the end of @self.
        It is assumed that other is an instant of Matrix
        """
        assert isinstance(other,Matrix)
        self._mat += [r[:] for r in other._mat]
        self._rows += other.rows()

    # overriding functions of ConstBaseMatrix

    def __getitem__(self,key):
        i,j = key
        return self._mat[i][j]

    # overriding functions of BaseMatrix

    def __setitem__(self,key,value):
        i,j = key
        self._mat[i][j] = value

class MatrixClass(BaseMatClass):
    """
    A metaclass for the matrix type of a Matrix object.
    """

    def __init__(self, field):
        self._field = field

    # must implement

    def zero_mat(self, rows, cols):
        """
        Returns a rows-by-cols zero matrix of the given matrix class.
        """
        return Matrix(self._field, rows, cols)

    def rand(self, rows, cols):
        """
        Returns a random rows-by-cols matrix of the given matrix class.
        """
        return Matrix(self._field, rows, cols, \
                      mat=[[self._field.rand() for j in range(cols)] for i in range(rows)])

    def vec_class(self):
        """
        Returns the vector class with with the matrix class works (or works best).
        """
        return VectorClass(self._field)

    def __eq__(self, other):
        return self._field is other._field
