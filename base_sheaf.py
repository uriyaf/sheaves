
from block_mat import *
from simplicial_complex import DeltaSet, DeltaSetMorphism

class BaseSheaf(object):
    """
    Base class for sheaf of vector spaces. Do not use implement unless you really know what you're doing.
    """

    def __init__(self, simp_complex):
        self._sc = simp_complex
        self._d = [None for i in range(simp_complex.dim()+2)] # the last place stores d[-1]
        self._d_quasi_inverse = [None for i in range(simp_complex.dim()+2)] # the last place stores d_quasi_inverse[-1]
        self._coh_bases = [None for i in range(simp_complex.dim()+2)]
        self._coh_dual_bases = [None for i in range(simp_complex.dim()+2)]
        
    
    # must implement.

    def rest_codim1(self,face,i):
        """
        Returns an matrix representing the restriction map from {face - vertex i} to {face}.
        """
        pass

    def mat_class(self):
        pass

    def face_dim(self,face):
        """
        Returns the dimension of the sheaf at the given face.
        """
        pass

    def _get_blocks(self,i):
        """
        Returns a Blocks object representing blocks labelled by the i-dimensional faces of the
        underlying simplicial complex, and with length equal to the dimension of the sheaf on these faces.
        Should work for negative i-s and also for i-s greater than the dimension of the base simplicial complex.
        """
        pass

    # functions you get for free

    def cochain_class(self):
        """
        returns the class of cochains of the given sheaf.
        The default is Cochain.
        The result must inherit from Cochain.
        """
        return Cochain

    def base_complex(self):
        """
        Returns the simplicial complex on which the sheaf is defined.
        """
        return self._sc
    
    def d(self,i):
        """
        Returns the i-th coboundary map as a SparseBlockMat object.
        """
        if -1 <= i and i <= self._sc.dim():
            if self._d[i] is not None:
                return self._d[i]
        blocks = {}
        row_blocks = self._get_blocks(i+1)
        col_blocks = self._get_blocks(i)
        if -1 <= i and i < self._sc.dim():
            for f in self._sc.get_faces(i+1):
                for j in range(0,i+2):
                    f1 = self._sc.subface_codim1(f,j) # f1 = f[:j]+f[j+1:] # omit j-th coordinate
                    if (f,f1) in blocks.keys():
                        #print("double face map", f1, f)
                        M = blocks[f,f1]
                        if j % 2 == 0:
                            blocks[f,f1] = M + self.rest_codim1(f,j)
                        else:
                            blocks[f,f1] = M - self.rest_codim1(f,j)
                    else:
                        if j % 2 == 0:
                            blocks[f,f1] = self.rest_codim1(f,j)
                        else:
                            blocks[f,f1] = -self.rest_codim1(f,j)
        res = SparseBlockMat(blocks, row_blocks, col_blocks, self.mat_class())
        if -1 <= i and i <= self._sc.dim():
            self._d[i] = res
        return res

    def d_quasi_inverse(self, i):
        """
        Returns a quasi-inverse of the i-th coboundary map.
        The output is a matrix having the same matrix class as the sheaf.
        """
        if -1 <= i and i <= self._sc.dim():
            if self._d_quasi_inverse[i] is not None:
                return self._d_quasi_inverse[i]
        res = self.d(i).quasi_inverse()
        if -1 <= i and i <= self._sc.dim():
            self._d_quasi_inverse[i] = res
        return res

    def rest(self, face, remove_verts):
        """
        Returns an matrix object representing the restriction map from {face - remove_verts} to {face}.
        @remove_verts are the indices of the vertices to be removed. They are consecutively,
        e.g., [0,1], will result in removing the 0th vertex and then removing the 1st vertex from the result.
        """
        #remove_verts = sorted(remove_verts, reverse=True) # best not to do it
        if len(remove_verts) == 0:
            return self.mat_class().id_mat(self.face_dim(face))
        g = face
        mat = self.rest_codim1(g,remove_verts[0])
        g = self._sc.subface_codim1(g,remove_verts[0])
        for i in remove_verts[1:]:
            mat = mat * self.rest_codim1(g,i)  
            g = self._sc.subface_codim1(g,i)
        return mat

    def forget_d(self):
        """
        Forgets the differential objects.
        """
        self._d = [None for i in range(simp_complex.dim()+2)] # the last place stores d[-1]

    def forget_d_quasi_inverse(self):
        """
        Forgets the quasi_inverses of the differentials.
        """
        self._d_qausi_inverse = [None for i in range(simp_complex.dim()+2)] # the last place stores d_quasi_inverse[-1]

    def forget_coh_bases(self):
        """
        Forgets the cohmology bases.
        """
        self._coh_bases = [None for i in range(simp_complex.dim()+2)]

    def forget_coh_dual_bases(self):
        """
        Forgets the cohmology dual bases.
        """
        self._coh_dual_bases = [None for i in range(simp_complex.dim()+2)]

    def total_dim(self,i):
        """
        Returns the total dimension of the i-dimensional faces.
        """
        return sum(self.face_dim(f) for f in self.base_complex().get_faces(i))

    def Euler_characteristic(self):
        """
        Returns the Euler characteristic of the sheaf.
        """
        res = 0
        sign = -1
        for d in range(-1, self._sc.dim()+1):
            res += sign * self.total_dim(d)
            sign = -sign
        return res

    def validate(self):
        """
        Checks whether the sheaf conditions are satisfied.
        """
        if self._sc.dim()<1:
            return True
        for dim in range(1,self._sc.dim()+1):
            for f in self._sc.get_faces(dim):
                for i in range(dim+1):
                    for j in range(i):
                        # print(f,i,j)
                        fi = self._sc.subface_codim1(f,i) #fi = f[:i] + f[i+1:]
                        fj = self._sc.subface_codim1(f,j) #fj = f[:j] + f[j+1:]
                        M = self.rest_codim1(f,i) * self.rest_codim1(fi,j)
                        N = self.rest_codim1(f,j) * self.rest_codim1(fj,i-1)
                        if M != N:
                            self._validation_failure = (f,i,j)
                            return False
        return True

    # functions involving cochains

    def cohomology_basis(self,i):
        """
        Returns an array of i-cocycles representing a basis for the i-th cohomology.
        """
        if i < -1 or i > self._sc.dim():
            return []
        if self._coh_bases[i] is not None:
            return self._coh_bases[i]
        raw = cohomology_basis(self.d(i),self.d(i-1)).get_rows()
        cochain_class = self.cochain_class()
        res = [cochain_class(self,i,v) for v in raw]
        self._coh_bases[i] = res
        return res

    def cohomology_dual_basis(self,i):
        """
        Returns an array of i-cocycles representing a dual basis for the i-th cohomology.
        (The cochains should be regarded as cochains of the dual sheaf.)
        """
        if i < -1 or i > self._sc.dim():
            return []
        if self._coh_dual_bases[i] is not None:
            return self._coh_dual_bases[i]
        raw = cohomology_dual_basis(self.d(i),self.d(i-1)).get_rows()
        cochain_class = self.cochain_class()
        res = [cochain_class(self,i,v) for v in raw]
        self._coh_dual_bases[i] = res
        return res

    def rand_cohomology_basis(self,i):
        """
        Returns an arrary of i-cocycles which represent a randomly chosen basis for the i-th cohomology.
        """
        if i < -1 or i > self._sc.dim():
            return []
        raw = rand_cohomology_basis(self.d(i),self.d(i-1)).get_rows()
        cochain_class = self.cochain_class()
        res = [cochain_class(self,i,v) for v in raw]
        return res

    def zero_cochain(self,i):
        """
        Returns a zero i-cochain.
        """
        return self.cochain_class()(self, i)

    def rand_cochain(self,i):
        """
        Returns a random i-cochain.
        """
        return self.cochain_class()(self, i, self.mat_class().vec_class().rand(self.total_dim(i)))

    def rand_coboundary(self,i):
        """
        Returns a random i-coboundary.
        """
        raw = self.d(i-1).apply(self.rand_cochain(i-1))
        return self.cochain_class()(self,i,raw)

    def rand_cocycle(self,i):
        """
        Returns a random i-cocycle. This may be time-expensive if dim H^i is large.
        """
        raw = self.d(i-1).apply(self.rand_cochain(i-1))
        return self.cochain_class()(self,i,raw) + \
               rand_linear_combination(self.cohomology_basis(i), self.zero_cochain(i))

    def source_for_d(self, cochain):
        """
        Takes a i-cochain, which is assumed to be a coboundary,
        and returns a (i-1)-cochain mapping to @cochain under d(i-1).
        This operation may be time-expensive.
        """
        i = cochain.get_dim()
        raw = self.d_quasi_inverse(i-1).apply(cochain._vec)
        return self.cochain_class()(self, i-1, raw)

    # pullback, pushforward and related functions

    def pullback(self, covering):
        """
        Returns the pullback of the sheaf along the covering of simplicial coverings given.
        @coverings is a DeltaSetMorphism object. 
        """
        assert covering.target == self._sc
        Y = covering.source
        rest = {}
        for d in range(Y.dim()+1):
            for i in range(d+1):
                for f in Y.get_faces(d):
                    rest[f,i] = self.rest_codim1(covering(f),i)
        return GeneralSheaf(Y, rest, mat_class=self.mat_class(), \
                            empty_face_dim=self.face_dim(self._sc.get_empty_face()))

    def pushforward(self, covering):
        """
        Returns the pushforward of the sheaf along the covering of simplicial coverings given.
        @covering is a DeltaSetMorphism object. It is assumed to be a topological covering.
        The sheaf @self has to be strictly a sheaf (and not an augmented sheaf).
        """
        assert covering.source == self._sc
        fibers = {}
        target = covering.target
        source = covering.source
        for f in target.get_all_faces():
            fibers[f] = []
        for f in source.get_all_faces():
            fibers[covering(f)].append(f)
        fiber_blocks = {}
        for f in target.get_all_faces():
            fiber_blocks[f] = Blocks([(g, self.face_dim(g)) for g in fibers[f]])
        rest = {}
        mat_class = self.mat_class()
        for d in range(0, target.dim()+1):
            for i in range(d+1):
                for f in target.get_faces(d):
                    blocks = {}
                    f1 = target.subface_codim1(f,i)
                    for g in fibers[f]:
                        g1 = source.subface_codim1(g,i)
                        blocks[g,g1] = self.rest_codim1(g,i)
                    rest[f,i] = SparseBlockMat(blocks, fiber_blocks[f], fiber_blocks[f1], mat_class=mat_class)
        return GeneralSheaf(target, rest, mat_class, len(fiber_blocks[target.get_empty_face()]))

    def unit_map(self, covering):
        """
        @covering is a DeltaSetMorphism object with target equaling to the sheaf's base complex.
        Returns a SheafMorphism representing the natural map from 'self' to 'cov_*cov^*(self)', a.k.a. the unit map.
        """
        assert self._sc == covering.target
        F = self.pullback(covering).pushforward(covering)
        unit = {}
        mat_class = self.mat_class()
        for f in self._sc.get_all_faces():
            M = mat_class.zero_mat(F.face_dim(f), self.face_dim(f))
            for i in range(F.face_dim(f)):
                M.set_to_one(i, i % self.face_dim(f))
            unit[f] = M
        return SheafMorphism(self,unit,F)

    def counit_map(self, covering):
        """
        @coving is a DeltaSetMorphism object with target equaling to the sheaf's base complex.
        Returns a SheafMorphism representing the natural map from 'cov_*cov^*(self)' to 'self', a.k.a. the counit map.
        """
        assert self._sc == covering.target
        F = self.pullback(covering).pushforward(covering)
        counit = {}
        mat_class = self.mat_class()
        for f in self._sc.get_all_faces():
            M = mat_class.zero_mat(self.face_dim(f), F.face_dim(f))
            for i in range(F.face_dim(f)):
                M.set_to_one(i % self.face_dim(f),i) 
            counit[f] = M
        return SheafMorphism(F,counit,self)
    
    
class Cochain(BlockVec):
    """
    Repersents a cochain of a sheaf. Also functions as a block vector.
    """
    def __init__(self, sheaf, dim, vec=None):
        self._sheaf = sheaf
        self._dim = dim
        if vec is None:
            BlockVec.__init__(self, sheaf._get_blocks(dim), vec, sheaf.mat_class().vec_class())
        else:
            BlockVec.__init__(self, sheaf._get_blocks(dim), vec)

    # overriding functions of BlockVec

    def clone(self):
        return Cochain(self._sheaf, self._dim, self._vec)

    def scale(self, scalar):
        """
        Returns the vector obtained from scaling @self by the given scalar.
        """
        return Cochain(self._sheaf, self._dim, self._vec.scale(scalar))

    def __add__(self, other):
        if not isinstance(other, Cochain):
            raise "Cannot add Cochain and non-Cochain" # other is not a block vector, so block structure is lost
        if self._sheaf != other._sheaf:
            raise "Sheaves are not compatible."
        if self._dim != other._dim:
            raise "Dimensions are not compatible."
        return Cochain(self._sheaf, self._dim, self._vec + other._vec)

    def __sub__(self, other):
        if not isinstance(other, Cochain):
            raise "Cannot add Cochain and non-Cochain" # other is not a block vector, so block structure is lost
        if self._sheaf != other._sheaf:
            raise "Sheaves are not compatible."
        if self._dim != other._dim:
            raise "Dimensions are not compatible."
        return Cochain(self._sheaf, self._dim, self._vec - other._vec)

    def __neg__(self):
        return Cochain(self._sheaf, self._dim, -self._vec)

    def __repr__(self):
        res = "%d-cochain:" % self._dim
        max_print = 10
        for k in self._blocks.keys():
            max_print -= 1
            if max_print < 0:
                res += "\n..."
                break
            res += "\n" + str(k)+": "+str(self.get_block(k))
        return res

    # additional functions

    def to_vector(self):
        """
        Returns a vector representing the cochain.
        """
        return self._vec.clone()

    def coboundary(self):
        i = self._dim
        raw = self._sheaf.d(i).apply(self._vec)
        return Cochain(self._sheaf, i+1, raw)

    def get_dim(self):
        return self._dim

    def get_sheaf(self):
        return self._sheaf





class SheafMorphism(object):
    """
    Implements a morphism between two sheaves over the same simplicial complex (i.e. DeltaSet object).
    """
    def __init__(self, source_sheaf, morphism, target_sheaf):
        """
        @morphism is a dictionary mapping faces of the base simplicial complex to matrices.
        """
        assert source_sheaf.base_complex() == target_sheaf.base_complex()
        self.source = source_sheaf
        self.target = target_sheaf
        self._morphism = morphism

    def validate(self):
        sc = self.source.base_complex()
        for d in range(sc.dim()+1):
            for f in sc.get_faces(d):
                for i in range(d+1):
                    f1 = sc.subface_codim1(f,i)
                    L1 = self._morphism[f] * self.source.rest_codim1(f,i)
                    L2 = self.target.rest_codim1(f,i) * self._morphism[f1]
                    if L1._mat != L2._mat:
                        return False
        return True

    def ker(self):
        """
        Returns a SheafMorphism representing the kernel of @self.
        """
        rest = {}
        sc = self.source._sc
        ker = {}
        for f in sc.get_all_faces():
            ker[f] = self._morphism[f].ker()
        for d in range(sc.dim() + 1):
            for f in sc.get_faces(d):
                Q = ker[f].quasi_inverse()
                for i in range(d+1):
                    f1 = sc.subface_codim1(f,i)
                    rest[f,i] = Q * self.source.rest_codim1(f, i) * ker[f1]
        F = GeneralSheaf(sc, rest, mat_class=self.source.mat_class(), \
                         empty_face_dim=ker[sc.get_empty_face()].rows())
        return SheafMorphism(F, ker, self.source)
    
    def coker(self):
        """
        Returns a SheafMorphism representing the cokernel of @self.
        """
        rest = {}
        mat_class = self.target.mat_class
        sc = self.target._sc
        coker = {}
        for f in sc.get_all_faces():
            coker[f] = self._morphism[f].coker()
        for d in range(sc.dim() + 1):
            coker_left_inverse = {g:coker[g].quasi_inverse() for g in sc.get_faces(d-1)}
            for f in sc.get_faces(d):
                for i in range(d+1):
                    f1 = sc.subface_codim1(f,i)
                    rest[f,i] = coker[f] * self.target.rest_codim1(f, i) * coker_left_inverse[f1]
        F = GeneralSheaf(sc, rest, mat_class=self.target.mat_class(), \
                         empty_face_dim=coker[sc.get_empty_face()].rows())
        return SheafMorphism(self.target, coker, F)

    def apply(self, cochain):
        """
        Applies the sheaf morphism to the given cochain and returns the result.
        """
        source = self.source
        assert source == cochain.get_sheaf()
        d = cochain.get_dim()
        res = self.target.zero_cochain(d) 
        for f in source.base_complex().get_faces(d):
            res.set_block(f, self._morphism[f].apply(cochain.get_block(f)))
        return res

    def __call__(self, cochain):
        """
        Applies the sheaf morphism to the given cochain and returns the result.
        """
        return self.apply(cochain)

class GeneralSheaf(BaseSheaf):
    """
    Implements a general sheaf over a field.
    """

    def __init__(self, simp_complex, restriction_maps, mat_class=None, empty_face_dim=None):
        """
        @simplicial_complex is a DeltaSet object
        @restriction_maps is a dictionary taking pairs (f,j) consisting of an
            i-dimensional face and j in {0,...,i} to a matrix of the given matrix class.
        @mat_class is the relevant matrix class. If omitted, it will be deduced from
            the restriction maps.
        @empty_face_dim is the dimension of the sheaf at the zeroeth face (typically 0).
            If omitted, it will be read from the restriction maps.
        Compatibility conditions are not verified, but you can use the validate method to check that the restrction maps are compatible.
        """
        BaseSheaf.__init__(self, simp_complex)
        if mat_class is None:
            for r in restriction_maps.values():
                mat_class = r.mat_class()
                break
        self._mat_class = mat_class
        self._rest = restriction_maps
        if empty_face_dim is None:
            for f in X.get_faces(0):
                empty_face_dim = restriction_maps[f,0].cols()
                break
        self._empty_face_dim = empty_face_dim
        self._blocks = [Blocks([(f, restriction_maps[f, 0].rows()) for f in simp_complex.get_faces(i)]) for i in range(simp_complex.dim()+1)] + \
                       [Blocks([(simp_complex.get_empty_face(), empty_face_dim)])]
        
    # implementation of functions from BaseSheaf 

    def rest_codim1(self,face,i):
        """
        Returns an matrix representing the restriction map from {face - vertex i} to {face}.
        """
        return self._rest[face,i]

    def mat_class(self):
        return self._mat_class

    def face_dim(self,face):
        """
        Returns the dimension of the sheaf at the given face.
        """
        if face == self._sc.get_empty_face():
            return self._empty_face_dim
        else:
            return self._rest[face,0].rows()

    def _get_blocks(self,i):
        """
        Returns a Blocks object representing blocks labelled by the i-dimensional faces of the
        underlying simplicial complex, and with length equal to the dimension of the sheaf on these faces.
        Should work for negative i-s and also for i-s greater than the dimension of the base simplicial complex.
        """
        if i < -1 or i > self._sc.dim():
            return Blocks([])
        return self._blocks[i]
