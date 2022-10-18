from base_sheaf import *


class StructureSheaf(BaseSheaf):
    """
    Represents a structure sheaf over a field.
    """

    def __init__(self, simp_complex, mat_class):
        BaseSheaf.__init__(self, simp_complex)
        self._mat_class = mat_class
        self._blocks = [Blocks([(f,1) for f in simp_complex.get_faces(i)]) for i in range(simp_complex.dim()+1)] + \
                       [Blocks([(simp_complex.get_empty_face(),0)])]
        self._rest = mat_class.const_id_mat(1)
        self._zero_rest = mat_class.const_zero_mat(1,0)

    # must implement for BaseSheaf

    def rest_codim1(self,face,i):
        """
        Returns an matrix representing the restriction map from {face - vertex i} to {face}.
        """
        if self._sc.face_dim(face) == 0:
            return self._zero_rest
        return self._rest 

    def mat_class(self):
        return self._mat_class

    def face_dim(self,face):
        """
        Returns the dimension of the sheaf at the given face.
        """
        if face != self._sc.get_empty_face():
            return 1
        return 0

    def _get_blocks(self,i):
        """
        Returns a Blocks object representing blocks labelled by the i-dimensional faces of the
        underlying simplicial complex, and with length equal to the dimension of the sheaf on these faces.
        """
        if i < -1 or i > self._sc.dim():
            return Blocks([])
        return self._blocks[i]

    # overriding functions from BaseSheaf

    def pullback(self, covering):
        """
        Returns the pullback of the sheaf along the covering of simplicial coverings given.
        @coverings is a DeltaSetMorphism object. 
        """
        assert covering.target == self._sc
        return StructureSheaf(covering.source, self._mat_class)



class ConstantSheaf(BaseSheaf):
    """
    Represents a constant sheaf over the field.
    """

    def __init__(self, simp_complex, dim, mat_class):
        BaseSheaf.__init__(self, simp_complex)
        self._mat_class = mat_class
        self._blocks = [Blocks([(f,dim) for f in simp_complex.get_faces(i)]) for i in range(simp_complex.dim()+1)] + \
                       [Blocks([(simp_complex.get_empty_face(),0)])]
        self._dim = dim
        self._rest = mat_class.const_id_mat(dim)
        self._zero_rest = mat_class.const_zero_mat(dim,0)

    # must implement for BaseSheaf

    def rest_codim1(self,face,i):
        """
        Returns an matrix representing the restriction map from {face - vertex i} to {face}.
        """
        if self._sc.face_dim(face) == 0:
            return self._zero_rest
        return self._rest

    def mat_class(self):
        return self._mat_class

    def face_dim(self,face):
        """
        Returns the dimension of the sheaf at the given face.
        """
        if face != self._sc.get_empty_face():
            return self._dim
        return 0

    def _get_blocks(self,i):
        """
        Returns a Blocks object representing blocks labelled by the i-dimensional faces of the
        underlying simplicial complex, and with length equal to the dimension of the sheaf on these faces.
        """
        if i < -1 or i > self._sc.dim():
            return Blocks([])
        return self._blocks[i]

    # overriding functions from BaseSheaf

    def pullback(self, covering):
        """
        Returns the pullback of the sheaf along the covering of simplicial coverings given.
        @coverings is a DeltaSetMorphism object. 
        """
        assert covering.target == self._sc
        return ConstantSheaf(covering.source, self._dim, self._mat_class)


class DirectSumSheaf(BaseSheaf):
    """
    A sheaf object represeting a direct sum of sheaves.
    The summand structure is exploited in various heavy computations to save time and memory.
    """
    def __init__(self, summands, simp_complex=None, mat_class=None):
        """
        @summands is a list of sheaves over the given simplicial complex having given matrix class.
        If omitted, @simp_complex and @mat_class are deduced from the list @summands (which must be nonempty).
        """
        if simp_complex is None:
            simp_complex = summands[0].base_complex()
        if mat_class is None:
            mat_class = summands[0].mat_class()
        for S in summands:
            assert S.base_complex() == simp_complex
            assert S.mat_class() == mat_class
        BaseSheaf.__init__(self, simp_complex)
        self._summands = summands[:]
        self._mat_class = mat_class
        subblocks = {}
        for f in simp_complex.get_all_faces():
            subblocks[f] = Blocks([(j, S.face_dim(f)) for j,S in enumerate(summands)])
        self._subblocks = subblocks
        self._blocks = [None]*(simp_complex.dim()+2) # the last index is the blocks in dim. -1.
        for d in range(-1, simp_complex.dim()+1):
            self._blocks[d] = Blocks([(f, len(subblocks[f])) for f in simp_complex.get_faces(d)])

    # must implement for BaseSheaf

    def rest_codim1(self,face,i):
        """
        Returns an matrix representing the restriction map from {face - vertex i} to {face}.
        """
        f1 = self._sc.subface_codim1(face,i)
        return SparseBlockMat({(j,j):S.rest_codim1(face,i) for j,S in enumerate(self._summands)}, \
                              self._subblocks[face], \
                              self._subblocks[f1], self._mat_class)

    def mat_class(self):
        return self._mat_class

    def face_dim(self,face):
        """
        Returns the dimension of the sheaf at the given face.
        """
        return len(self._subblocks[face])

    def _get_blocks(self,i):
        """
        Returns a Blocks object representing blocks labelled by the i-dimensional faces of the
        underlying simplicial complex, and with length equal to the dimension of the sheaf on these faces.
        Should work for negative i-s and also for i-s greater than the dimension of the base simplicial complex.
        """
        if i < -1 or i > self._sc.dim():
            return Blocks([])
        return self._blocks[i]

    # overriding functions from BaseSheaf

    def d_quasi_inverse(self, i):
        """
        Returns a quasi-inverse of the i-th coboundary map.
        The output is a matrix having the same matrix class as the sheaf.
        """
        sc = self._sc
        if -1 <= i and i <= sc.dim():
            if self._d_quasi_inverse[i] is not None:
                return self._d_quasi_inverse[i]
        Q = [S.d_quasi_inverse(i) for S in self._summands]
        res_blocks = {}
        for f in sc.get_faces(i+1):
            for g in sc.get_faces(i):
                mat_blocks = {}
                for j,S in enumerate(self._summands):
                    a,b = S._get_blocks(i+1).interval(f)
                    c,d = S._get_blocks(i).interval(g)
                    mat_blocks[j,j] = Q[j].get_rect(c, d, a, b)
                res_blocks[g,f] = SparseBlockMat(mat_blocks, self._subblocks[g], self._subblocks[f], self._mat_class)
        res = SparseBlockMat(res_blocks, self._get_blocks(i), self._get_blocks(i+1), self._mat_class)
        if -1 <= i and i <= self._sc.dim():
            self._d_quasi_inverse[i] = res
        return res

    def cohomology_basis(self,i):
        """
        Returns an array of i-cocycles representing a basis for the i-th cohomology.
        """
        if i < -1 or i > self._sc.dim():
            return []
        if self._coh_bases[i] is not None:
            return self._coh_bases[i]
        bases = [S.cohomology_basis(i) for S in self._summands]
        res = []
        for j,B in enumerate(bases):
            base_vectors = [Cochain(self,i) for k in range(len(B))]
            for c,b in zip(base_vectors, B):
                for f in self._sc.get_faces(i):
                    of1 = self._blocks[i].interval(f)[0]
                    of2 = self._subblocks[f].interval(j)[0]
                    c.set_segment(of1 + of2, b(f))
            res += base_vectors
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
        bases = [S.cohomology_dual_basis(i) for S in self._summands]
        res = []
        for j,B in enumerate(bases):
            base_vectors = [Cochain(self,i) for k in range(len(B))]
            for c,b in zip(base_vectors, B):
                for f in self._sc.get_faces(i):
                    of1 = self._blocks[i].interval(f)[0]
                    of2 = self._subblocks[f].interval(j)[0]
                    c.set_segment(of1 + of2, b(f))
            res += base_vectors
        self._coh_dual_bases[i] = res
        return res

    def rand_cohomology_basis(self,i):
        """
        Returns an arrary of i-cocycles which represent a randomly chosen basis for the i-th cohomology.
        """
        if i < -1 or i > self._sc.dim():
            return []
        B = self.cohomology_basis(i)
        M = self._mat_class.rand_inv_mat(len(B))
        zero = Cochain(self, i)
        return [linear_combination(r, B, zero) + self.rand_coboundary(i) for r in M.get_rows()]

    def source_for_d(self, cochain):
        """
        Takes a i-cochain, which is assumed to be a coboundary,
        and returns a (i-1)-cochain mapping to @cochain under d(i-1).
        This operation may be time-expensive.
        """
        # this is typically faster than applying d_quasi_inverse directly
        return self.merge_cochains([S.source_for_d(c) for S,c in zip(self._summands, \
                                                                     self.decompose_cochain(cochain))])

    # additional methods

    def decompose_cochain(self, cochain):
        """
        Returns an array of cochains forming the components of the given
        cochain relative to the direct sum structure
        """
        assert self == cochain.get_sheaf()
        d = cochain.get_dim()
        res = [S.zero_cochain(d) for S in self._summands]
        for f in self._sc.get_faces(d):
            of = self._get_blocks(d).interval(f)[0]
            subblock = self._subblocks[f]
            for j, S in enumerate(self._summands):
                a, b = subblock.interval(j)
                res[j].set_block(f, cochain.get_segment(of + a, of + b))
        return res

    def merge_cochains(self, cochains, dim=None):
        """
        Takes a list of cochains such that the i-th cochain is a cochain of the i-th summand
        and returns a cochain of @self whose components are the given cochains.
        @dim is the common dimension of the given cochains. If omitted, it is read from the first cochain.
        """
        if dim is None:
            dim = cochains[0].get_dim()
        res = Cochain(self, dim)
        for f in self._sc.get_faces(dim):
            of1 = self._get_blocks(dim).interval(f)[0]
            subblock = self._subblocks[f]
            for j, S in enumerate(self._summands):
                of2 = subblock.interval(j)[0]
                res.set_segment(of1 + of2, cochains[j](f))
        return res
        


def subsheaf_generated_by_1cochains(sheaf, cochains):
    """
    Takes a sheaf and an array of 1-cochains and returns a SheafMorphism object
    representing an injective sheaf morphism with target @sheaf such that its
    image is the subsheaf generated by the given 1-cochains.
    """
    for c in cochains:
        assert c.get_dim() == 1
        assert c.get_sheaf() == sheaf
    sc = sheaf.base_complex()
    mat_class = sheaf.mat_class()
    # constructing dictionary underlying the resulting morphism
    inc = {} 
    for d in range(-1, sc.dim()+1):
        for f in sc.get_faces(d):
            M = mat_class.zero_mat( len(cochains) * d * (d+1) // 2, sheaf.face_dim(f) )
            r = 0
            for c in cochains:
                for i in range(d+1):
                    for j in range(i):
                        verts_to_remove = [v for v in range(0,j)] + [v for v in range(j+1,i)] + [v for v in range(i+1,d+1)]
                        edge = sc.subface(f, verts_to_remove[::-1])
                        M.set_row(r, sheaf.rest(f, verts_to_remove[::-1]).apply(c(edge)) )
                        r += 1
            inc[f] = M.Gauss_elim().transpose()
    # constructing the restriction maps
    rest = {}
    for d in range(0, sc.dim()+1):
        for f in sc.get_faces(d):
            Q = inc[f].quasi_inverse() # left inverse is more appropriate
            for i in range(d+1):
                f1 = sc.subface_codim1(f,i)
                rest[f,i] = Q * sheaf.rest_codim1(f,i) * inc[f1]
    # finish
    C = GeneralSheaf(sc, rest, mat_class, 0)
    return SheafMorphism(C, inc, sheaf)
