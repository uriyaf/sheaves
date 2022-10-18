
# contains test functions for the various classes in the project
# each funtions performs a bunch of tests, and returns a boolean indicating of all the tests were successful and
# a list of the failed tests.
# Many tests are randomized, and there is a very small chance that they fail.
# The tests also feature how one might use the different classes of the project.

from base_mat import *
from bit_mat import *
from block_mat import *
from matrix import *
from fields import *
from base_sheaf import *
from simplicial_complex import *
from sheaf_tools import *
from modification_process import *
from ramanujan_complex import *
from random import random


def _test_DirectSumSheaf(mat_class):
    res = []
    #
    cov = torus_2dim_morphism(grid=(3,3), covering_factors=(3,3))
    X = cov.source
    Y = cov.target
    #
    F = StructureSheaf(Y, mat_class)
    S = DirectSumSheaf([F,F,F])
    G = ConstantSheaf(Y, 3, mat_class)
    res.append( (0, S.validate() ) )
    bases = [S.cohomology_basis(i) for i in range(-1,4)]
    res.append( (1, [len(B) for B in bases] == [0,3,6,3,0]) )
    dual_bases = [S.cohomology_dual_basis(i) for i in range(-1,4)]
    res.append( (2, [len(D) for D in dual_bases] == [0,3,6,3,0]) )
    flag = True
    for B,D in zip(bases, dual_bases):
        flag = flag and (len(B)==len(D))
        for b in B:
            flag = flag and b.coboundary().is_zero()
        M = mat_class.zero_mat(len(B), len(D))
        for i,b in enumerate(B):
            for j,d in enumerate(D):
                M[i,j] = b*d
        #print(M)
        flag = flag and (M.inverse() is not None)
    res.append( (2.5, flag) )
    res.append( (3, S.mat_class() == mat_class) )
    res.append( (4, S.base_complex() == Y) )
    res.append( (4.1, S.d(1) == G.d(1)) )
    c = S.rand_coboundary(1)
    s = S.source_for_d(c)
    res.append( (4.2, s.coboundary() == c) ) 
    #
    G = StructureSheaf(X, mat_class).pushforward(cov)
    S = DirectSumSheaf([F,G], simp_complex=Y, mat_class=mat_class)
    res.append( (5, S.validate() ) )
    bases = [S.rand_cohomology_basis(i) for i in range(-1,4)]
    res.append( (6, [len(B) for B in bases] == [0,2,4,2,0]) )
    dual_bases = [S.cohomology_dual_basis(i) for i in range(-1,4)]
    res.append( (7, [len(D) for D in dual_bases] == [0,2,4,2,0]) )
    flag = True
    for B,D in zip(bases, dual_bases):
        flag = flag and (len(B)==len(D))
        for b in B:
            flag = flag and b.coboundary().is_zero()
        M = mat_class.zero_mat(len(B), len(D))
        for i,b in enumerate(B):
            for j,d in enumerate(D):
                M[i,j] = b*d
        #print(M)
        flag = flag and (M.inverse() is not None)
    res.append( (7.5, flag) )
    res.append( (8, S.mat_class() == mat_class) )
    res.append( (9, S.base_complex() == Y) )
    #
    c = S.rand_cochain(1)
    L = S.decompose_cochain(c)
    res.append( (9.1, len(L) == 2) )
    c2 = S.merge_cochains(L,1)
    res.append( (9.2, c == c2) )
    #
    S = DirectSumSheaf([], simp_complex=Y, mat_class=mat_class)
    res.append( (10, S.validate() ) )
    res.append( (11, [len(S.rand_cohomology_basis(i)) for i in range(-1,4)] == [0,0,0,0,0]) )
    res.append( (12, [len(S.cohomology_dual_basis(i)) for i in range(-1,4)] == [0,0,0,0,0]) )
    res.append( (13, S.mat_class() == mat_class) )
    res.append( (14, S.base_complex() == Y) )
    res.append( (15, [S.total_dim(i) for i in range(-1,4)] == [0,0,0,0,0]) )
    #
    return (False not in [y for x,y in res]), [x for x,y in res if y == False] 
    

def _test_SheafMorphism(mat_class):
    # The methods ker and coker of SheafMorphism are already tested in
    # _test_sheaf_unit_counit_ker_and_coker, _test_push_forward_and_pullback
    # Here we mostly test the apply function.
    res = []
    cov = torus_2dim_morphism(grid=(3,3), covering_factors=(3,3))
    X = cov.source
    Y = cov.target
    F = StructureSheaf(Y, mat_class)
    unit = F.unit_map(cov)
    res.append( (1, unit.validate()) )
    a = unit.source.rand_cochain(1)
    b = unit.source.rand_cochain(1)
    res.append( (2, not a.is_zero()) ) # make sure that the test is nontrivial
    res.append( (3, unit(a).coboundary() == unit(a.coboundary())) )
    res.append( (4, unit(a+b) == unit(a)+unit(b)) )
    res.append( (5, isinstance(unit(a), Cochain) ) )
    return (False not in [y for x,y in res]), [x for x,y in res if y == False] 


def _test_sheaf_unit_counit_ker_and_coker():
    res = []
    #
    F3 = PrimeField(3) 
    cov = torus_2dim_morphism(grid=(3,3), covering_factors=(3,3))
    X = cov.source
    Y = cov.target
    #
    F = StructureSheaf(Y, MatrixClass(F3))
    unit = F.unit_map(cov)
    res.append( (1, unit.validate() ) )
    res.append( (1.1, unit.source is F) )
    res.append( (1.2, unit.target.validate() ) )
    ker = unit.ker()
    res.append( (1.5, unit.source is ker.target) )
    res.append( (1.6, ker.source.validate()) )
    res.append( (1.7, [ker.source.total_dim(i) for i in range(3)] == [0,0,0]) )
    coker = unit.coker()
    res.append( (1.8, unit.target is coker.source) )
    res.append( (1.81, coker.target.validate()) )
    res.append( (1.82, coker.validate()) )
    res.append( (1.83, coker.target.Euler_characteristic() == 0) )
    counit = F.counit_map(cov)
    res.append( (2, counit.validate() ) )
    res.append( (2.1, counit.target is F) )
    res.append( (2.2, counit.source.validate() ) )
    ker = counit.ker()
    res.append( (2.5, counit.source is ker.target) )
    res.append( (2.6, ker.source.validate()) )
    res.append( (2.7, ker.validate() ) )
    res.append( (2.8, ker.source.Euler_characteristic() == 0) )
    #
    F = ConstantSheaf(Y, 2, BitMatClass(60))
    unit = F.unit_map(cov)
    res.append( (3, unit.validate() ) )
    res.append( (3.1, unit.source is F) )
    res.append( (3.2, unit.target.validate() ) )
    coker = unit.coker()
    res.append( (3.8, unit.target is coker.source) )
    res.append( (3.81, coker.target.validate()) )
    res.append( (3.82, coker.validate()) )
    res.append( (3.83, [len(coker.target.rand_cohomology_basis(i)) for i in range(-1,4)]==[0,0,0,0,0]) )
    counit = F.counit_map(cov)
    res.append( (4, counit.validate() ) )
    res.append( (4.1, counit.target is F) )
    res.append( (4.2, counit.source.validate() ) )
    ker = counit.ker()
    res.append( (5, counit.source is ker.target) )
    res.append( (6, ker.source.validate()) )
    res.append( (7, ker.validate() ) ) 
    res.append( (8, [len(ker.source.cohomology_basis(i)) for i in range(3)]==[0,0,0]) )
    #
    return (False not in [y for x,y in res]), [x for x,y in res if y == False] 


def _test_pushforward_and_pullback():
    res = []
    #
    F3 = PrimeField(3) 
    cov = torus_2dim_morphism(grid=(3,3), covering_factors=(3,3))
    X = cov.source
    Y = cov.target
    #
    F = StructureSheaf(Y, MatrixClass(F3))
    G = F.pullback(cov)
    res.append( (0, G.validate()) )
    res.append( (1, [len(G.cohomology_basis(i)) for i in range(3)] == [1,2,1]) )
    print(res[-1])
    #
    H = G.pushforward(cov)
    res.append( (2, H.validate()) )
    res.append( (3, [len(G.cohomology_basis(i)) for i in range(4)] == [1,2,1,0]) )
    print(res[-1])
    #
    B = BitMatClass(240)
    cov = torus_ndim_morphism(grid=(2,2,1), factors=(2,2,2))
    X = cov.source
    Y = cov.target
    #
    F = ConstantSheaf(Y, 2, B)
    G = F.pullback(cov)
    res.append( (4, G.validate()) )
    res.append( (5, [len(G.cohomology_basis(i)) for i in range(-1,4)] == [0,2,6,6,2]) )
    print(res[-1])
    #
    G = ConstantSheaf(X, 2, B)
    H = G.pushforward(cov)
    res.append( (6, H.validate()) )
    res.append( (7, [len(G.cohomology_basis(i)) for i in range(5)] == [2,6,6,2,0]) )
    print(res[-1])
    #
    return (False not in [y for x,y in res]), [x for x,y in res if y == False] 

    

def _test_subsheaf_generated_by_1cochains():
    res = []
    #
    F7 = PrimeField(7)
    X = torus_2dim(grid = (2,3))
    S = ConstantSheaf(X,4,MatrixClass(F7))
    c = S.rand_cochain(1)
    inc = subsheaf_generated_by_1cochains(S,[c,c])
    F = inc.source
    X = F.base_complex()
    L0 = [F.face_dim(f) for f in X.get_faces(0)]
    L1 = [F.face_dim(f) for f in X.get_faces(1)]
    L2 = [F.face_dim(f) for f in X.get_faces(2)]
    Lneg = [F.face_dim(f) for f in X.get_faces(-1)]
    res.append( (1, inc.target is S) )
    res.append( (2, max(L0) == 0 and min(L0) == 0) )
    res.append( (3, max(L1) == 1 and min(L1) == 1) )
    res.append( (4, max(L2) == 3 and min(L2) == 3) )
    res.append( (5, max(Lneg) == 0 and min(Lneg) == 0) )
    #
    X = torus_2dim(grid = (3,3))
    S = ConstantSheaf(X,10,BitMatClass(60))
    c1 = S.rand_cocycle(1)
    c2 = S.rand_cocycle(1)
    inc = subsheaf_generated_by_1cochains(S,[c1,c2,c1+c2])
    F = inc.source
    X = F.base_complex()
    L0 = [F.face_dim(f) for f in X.get_faces(0)]
    L1 = [F.face_dim(f) for f in X.get_faces(1)]
    L2 = [F.face_dim(f) for f in X.get_faces(2)]
    Lneg = [F.face_dim(f) for f in X.get_faces(-1)]
    res.append( (11, inc.target is S) )
    res.append( (12, max(L0) == 0 and min(L0) == 0) )
    res.append( (13, max(L1) == 2 and min(L1) == 2) )
    res.append( (14, max(L2) == 4 and min(L2) == 4) )
    res.append( (15, max(Lneg) == 0 and min(Lneg) == 0) )
    return (False not in [y for x,y in res]), [x for x,y in res if y == False] 

def _test_Cochain(mat_class):
    res = []
    X = torus_2dim(grid = (4,3))
    S = StructureSheaf(X, mat_class)
    # testing cochain arithmetic
    a = S.rand_cochain(0)
    b = S.rand_cochain(0)
    c = S.rand_cochain(0)
    res.append( (1, a != b) )
    res.append( (2, not a == b) )
    res.append( (2.5, a == a.clone()) )
    res.append( (3, a+b == b+a) )
    res.append( (4, (a-b) == -(b-a) ) )
    res.append( (5, isinstance(a+b, Cochain)) )
    res.append( (6, isinstance(a-b, Cochain)) )
    res.append( (7, isinstance(-a, Cochain)) )
    res.append( (7, isinstance(a.clone(), Cochain)) )
    res.append( (8, a*b + a*c == a*(b+c) ) ) # may fail for bit matrices
    print(a)
    res.append( (9, a.coboundary().coboundary().is_zero()) )
    res.append( (10, a.coboundary() + b.coboundary() == (a+b).coboundary()) )
    s = a.vec_class().rand(1)[0] # fishy way for getting a scalar
    res.append( (11, isinstance(a.scale(s), Cochain)) )
    res.append( (12, not a.scale(s).coboundary() != a.coboundary().scale(s)) )
    res.append( (13, not isinstance(a.to_vector(), Cochain) ) )
    d = Cochain(S,0,a.to_vector())
    res.append( (14, a==d) )
    for f in X.get_faces(0):
        a.set_block(f, b(f))
    res.append( (15, a==b) )
    # other extremeties
    a = Cochain(S,-1)
    res.append( (16, a.is_zero() ) )
    res.append( (17, a.coboundary().is_zero() ) )
    a = Cochain(S,4)
    res.append( (16, a.is_zero() ) )
    res.append( (17, a.coboundary().is_zero() ) )
    print(a)
    b = Cochain(S,4)
    res.append( (18, a+b == -((-a)-b)) )
    # random functions
    c = S.rand_cocycle(1)
    res.append( (19, not c.is_zero()) )
    res.append( (20, c.coboundary().is_zero()) )
    return (False not in [y for x,y in res]), [x for x,y in res if y == False] 

def _test_ConstantSheaf(mat_class):
    res = []
    # basic cohomology computations for 2-dim torus
    X = torus_2dim(grid = (3,3))
    S = ConstantSheaf(X, 3, mat_class)
    res.append( (0, S.validate()) )
    res.append( (1, len(S.cohomology_basis(0))==3) )
    B1 = S.cohomology_basis(1)
    res.append( (2, len(B1)==6) )
    res.append( (3, len(S.cohomology_basis(2))==3) )
    res.append( (4, len(S.cohomology_basis(3))==0) )
    res.append( (5, len(S.cohomology_basis(-1))==0) )
    res.append( (6, S.cohomology_basis(1) == B1) )
    res.append( (7, S.rand_cohomology_basis(1) != B1) )
    res.append( (8, len(S.rand_cohomology_basis(1)) == 6) )
    res.append( (9, len(S.cohomology_dual_basis(0)) == 3) )
    res.append( (10, len(S.cohomology_dual_basis(1)) == 6) )
    res.append( (11, len(S.cohomology_dual_basis(2)) == 3) )
    res.append( (12, S.d(1).apply(B1[0]).is_zero()) )
    f = X.rand_face(1)
    f1 = X.subface_codim1(f,0)
    res.append( (13, S.rest(f,[0,0]) == S.rest_codim1(f,0) * S.rest_codim1(f1,0)) )
    f = X.rand_face(2)
    f1 = X.subface_codim1(f,1)
    res.append( (14, S.rest(f,[1,0]) == S.rest_codim1(f,1) * S.rest_codim1(f1,0)) )
    # some cohomology computations for 3-dim torus
    X = torus_ndim(grid = (2,1,2))
    S = ConstantSheaf(X, 3, mat_class)
    res.append( (21, len(S.cohomology_basis(0))==3) )
    B1dual = S.cohomology_dual_basis(1)
    res.append( (22, len(B1dual)==9) )
    res.append( (23, len(S.cohomology_basis(2))==9) )
    res.append( (24, len(S.cohomology_basis(3))==3) )
    res.append( (25, len(S.cohomology_basis(-1))==0) )
    res.append( (26, S.cohomology_dual_basis(1) == B1dual) )
    res.append( (27, len(S.cohomology_basis(1)) == 9) )
    res.append( (28, len(S.rand_cohomology_basis(1)) == 9) )
    res.append( (29, len(S.cohomology_dual_basis(0)) == 3) )
    res.append( (30, len(S.cohomology_dual_basis(1)) == 9) )
    res.append( (31, len(S.cohomology_dual_basis(3)) == 3) )
    #
    B1 = S.rand_cohomology_basis(1)
    M = mat_class.zero_mat(len(B1), len(B1dual))
    for i,c in enumerate(B1):
        for j, d in enumerate(B1dual):
            M[i,j] = c*d
    res.append( (32, M.inverse() is not None) )
    #
    b = S.rand_cocycle(1)
    res.append( (32.5, b.coboundary().is_zero()) )
    V = mat_class.vec_class().zero_vec(len(B1dual))
    for j,d in enumerate(B1dual):
        V[j] = b*d
    res.append( (33, not V.is_zero()) )
    res.append( (34, not b.is_zero() ) )
    #
    z = S.zero_cochain(2)
    res.append( (35, z.is_zero()) )
    #
    b = S.rand_coboundary(2)
    s = S.source_for_d(b)
    res.append( (36, s.coboundary() == b) )
    s2 = S.source_for_d(b)
    res.append( (37, s == s2) )
    #
    return (False not in [y for x,y in res]), [x for x,y in res if y == False] 

def _test_StructureSheaf(mat_class):
    res = []
    # basic cohomology computations for 2-dim torus
    X = torus_2dim(grid = (3,3))
    S = StructureSheaf(X, mat_class)
    res.append( (0, S.validate()) )
    res.append( (1, len(S.cohomology_basis(0))==1) )
    B1 = S.cohomology_basis(1)
    res.append( (2, len(B1)==2) )
    res.append( (3, len(S.cohomology_basis(2))==1) )
    res.append( (4, len(S.cohomology_basis(3))==0) )
    res.append( (5, len(S.cohomology_basis(-1))==0) )
    res.append( (6, S.cohomology_basis(1) == B1) )
    res.append( (7, S.rand_cohomology_basis(1) != B1) )
    res.append( (8, len(S.rand_cohomology_basis(1)) == 2) )
    res.append( (9, len(S.cohomology_dual_basis(0)) == 1) )
    res.append( (10, len(S.cohomology_dual_basis(1)) == 2) )
    res.append( (11, len(S.cohomology_dual_basis(2)) == 1) )
    res.append( (12, S.d(1).apply(B1[0]).is_zero()) )
    # some cohomology computations for 3-dim torus
    X = torus_ndim(grid = (1,1,1))
    S = StructureSheaf(X, mat_class)
    res.append( (21, len(S.cohomology_basis(0))==1) )
    B1dual = S.cohomology_dual_basis(1)
    res.append( (22, len(B1dual)==3) )
    res.append( (23, len(S.cohomology_basis(2))==3) )
    res.append( (24, len(S.cohomology_basis(3))==1) )
    res.append( (25, len(S.cohomology_basis(-1))==0) )
    res.append( (26, S.cohomology_dual_basis(1) == B1dual) )
    res.append( (27, len(S.cohomology_basis(1)) == 3) )
    res.append( (28, len(S.rand_cohomology_basis(1)) == 3) )
    res.append( (29, len(S.cohomology_dual_basis(0)) == 1) )
    res.append( (30, len(S.cohomology_dual_basis(1)) == 3) )
    res.append( (31, len(S.cohomology_dual_basis(3)) == 1) )
    #
    B1 = S.rand_cohomology_basis(1)
    M = mat_class.zero_mat(len(B1), len(B1dual))
    for i,c in enumerate(B1):
        for j, d in enumerate(B1dual):
            M[i,j] = c*d
    res.append( (32, M.inverse() is not None) )
    # back to 2-dim torus
    X = torus_2dim(grid = (3,3))
    S = StructureSheaf(X, mat_class)
    B1 = S.rand_cohomology_basis(1)
    B1dual = S.cohomology_dual_basis(1)
    b = S.rand_coboundary(1)
    V = mat_class.vec_class().zero_vec(len(B1dual))
    for j,d in enumerate(B1dual):
        V[j] = b*d
    res.append( (33, V.is_zero()) )
    res.append( (34, not b.is_zero() ) )
    #
    z = S.zero_cochain(2)
    res.append( (35, z.is_zero()) )
    return (False not in [y for x,y in res]), [x for x,y in res if y == False] 


def _rand_SparseBlockMat(mat_class, row_blocks, col_blocks, p = 0.5):
    blocks = {}
    for rk in row_blocks.keys():
        for ck in col_blocks.keys():
            if random() < p:
                blocks[rk,ck] = mat_class.rand(row_blocks.block_len(rk), col_blocks.block_len(ck))
    return SparseBlockMat(blocks, row_blocks, col_blocks, mat_class)

def _test_SparseBlockMat(mat_class, blocks1 = None, blocks2 = None, blocks3 = None):
    if blocks1 is None:
        blocks1 = Blocks([('a',13), ('b',11), ('c', 15)])
    if blocks2 is None:
        blocks2 = Blocks([('a',2), ('e',16), ('k', 9)])
    if blocks3 is None:
        blocks3 = Blocks([(3,7), (9,0), ('x',14), ('y', 14)])
    blocks2_prime = Blocks([(1, len(blocks2))])
    res = []
    # test arithmetic operations
    a = _rand_SparseBlockMat(mat_class, blocks1, blocks2)
    b = _rand_SparseBlockMat(mat_class, blocks1, blocks2)
    c = _rand_SparseBlockMat(mat_class, blocks2, blocks3)
    a_prime = _rand_SparseBlockMat(mat_class, blocks1, blocks2_prime)
    c_prime = _rand_SparseBlockMat(mat_class, blocks2_prime, blocks3)
    A = a.forget_blocks()
    B = b.forget_blocks()
    C = c.forget_blocks()
    res.append( (1, a+b == A+B) )
    res.append( (1.5, a+B == A+b) )
    res.append( (2, a-b == A-B) )
    res.append( (2.5, a-B == A-b) )
    res.append( (3, a*c == A*C) )
    res.append( (3.5, A*c == a*C) )
    res.append( (4, -a == -A) )
    res.append( (5, a+b == b+a) )
    res.append( (6, (a-b)*c == a*c - b*c) )
    res.append( (7, a+b == b+a) )
    res.append( (8, a-b == -(b-a)) )
    res.append( (9, a_prime+b == b+a_prime) ) # blocks do not match
    res.append( (10, a_prime-b == -(b-a_prime)) ) # blocks do not match
    res.append( (11, (a-b)*c_prime == a*c_prime - b*c_prime) ) # blocks do not match
    # test eq and neq
    a = _rand_SparseBlockMat(mat_class, blocks1, blocks2)
    b = _rand_SparseBlockMat(mat_class, blocks1, blocks2)
    c = _rand_SparseBlockMat(mat_class, blocks1, blocks2_prime)
    d = _rand_SparseBlockMat(mat_class, blocks2, blocks3)
    A = a.forget_blocks()
    res.append( (12, a==a) )
    res.append( (12.1, a==A) )
    res.append( (12.2, A==a) )
    res.append( (12.3, not a!=A) )
    res.append( (12.4, not A!=a) )
    res.append( (12.5, a==a.clone()) )
    res.append( (12.7, a==a.to_mutable()) )
    res.append( (13, not a!=a) )
    res.append( (14, not a==b) )
    res.append( (15, a!=b) )
    res.append( (16, not a==c) )
    res.append( (17, a!=c) )
    res.append( (18, not a==d) )
    res.append( (19, a!=d) )
    res.append( (20, a == a.forget_blocks()) )
    res.append( (21, not a != a.forget_blocks()) )
    res.append( (22, a.forget_blocks() == a) )  
    res.append( (23, not a.forget_blocks() != a) )
    res.append( (24, (a-a).is_zero()) )
    res.append( (25, not a.is_zero()) )
    # test Gauss elimination
    a = _rand_SparseBlockMat(mat_class, blocks1, blocks2)
    e = _rand_SparseBlockMat(mat_class, blocks1, blocks3)
    g,l,e_new = a.Gauss_elim_extended_extras(e)
    A = a.forget_blocks()
    E = e.forget_blocks()
    g,l,e_new = a.Gauss_elim_extended_extras(e)
    G,L,E_new = A.Gauss_elim_extended_extras(E)
    res.append( (26, g==G) )
    res.append( (27, l==L) )
    res.append( (28, e_new==E_new) )
    g,l = a.Gauss_elim_extended()
    G,L = A.Gauss_elim_extended()
    res.append( (29, g==G) )
    res.append( (30, l==L) )
    g = a.Gauss_elim()
    G = A.Gauss_elim()
    res.append( (31, g==G) )
    # test get functions
    a = _rand_SparseBlockMat(mat_class, blocks1, blocks2)
    r = randint(0,a.rows()-1)
    c = randint(0,a.cols()-1)
    A = a.forget_blocks()
    res.append( (32, a.get_row(r) == A.get_row(r)) )
    res.append( (33, a.get_col(c) == A.get_col(c)) )
    res.append( (34, a.get_entry(r,c) == A[r,c]) )
    # test apply
    a = _rand_SparseBlockMat(mat_class, blocks1, blocks2)
    v = mat_class.vec_class().rand(len(blocks2))
    u = BlockVec(blocks2, v)
    res.append( (35, a.apply(v) == a.forget_blocks().apply(v)) )
    res.append( (36, a.apply(u) == a.apply(v)) )
    res.append( (37, not isinstance(a.apply(u), BlockVec) ) )
    res.append( (38, a.apply(v).vec_class() == mat_class.vec_class()) )
    res.append( (39, a.apply_and_return_as_BlockVec(v) == a.forget_blocks().apply(v)) )
    res.append( (40, isinstance(a.apply_and_return_as_BlockVec(v), BlockVec) ) )
    # test matrix class behavior
    a = _rand_SparseBlockMat(mat_class, blocks1, blocks2)
    res.append( (100, a.mat_class() == a.clone().mat_class()) )
    res.append( (101, a.mat_class() == a.to_mutable().mat_class()) )
    res.append( (101, a.mat_class() == a.forget_blocks().mat_class()) )
    res.append( (102, a.validate()) )
    #
    return (False not in [y for x,y in res]), [x for x,y in res if y == False] 


def _test_BlockVec_visually(vec_class, blocks1=None, blocks2=None):
    """
    It is assumed that len(blocks1) == len(blocks2) > 10 if they are not None.
    """
    if blocks1 is None:
        blocks1 = Blocks([('a',4), ('b',0), ('c',9)])
    if blocks2 is None:
        blocks2 = Blocks(t for t in [('x',5), ('b',2), ('c',2), ('y',4)])
    # testing with matching blocks
    print("=============================")
    a = BlockVec(blocks1, vec_class.rand(len(blocks1)))
    b = BlockVec(blocks1, vec_class.rand(len(blocks1)))
    z = BlockVec(blocks1, vec_class = vec_class)
    print("a = ", a)
    print("b = ", b)
    print("z = ", z)
    print("len(a) = ", len(a))
    print("len(z) = ", len(z))
    print("a.clone() = ", a.clone())
    print("a[1] = ", a[1], "=", a.get_entry(1))
    print("Setting a[1] = b[1] and a[2] = b[2]")
    a.set_entry(1,b.get_entry(1))
    a[2] = b[2]
    print("a = ", a)
    print("b = ", b)
    print("b[3:5] = ", b.get_segment(3,5))
    print("setting a[5:7] to b[3:5]")
    a.set_segment(5, b.get_segment(3,5))
    print("a =   ", a)
    print("b =   ", b)
    print("a+b = ", a+b)
    print("a-b = ", a-b)
    print("-a  = ", a)
    print("a*b = ", a*b)
    print("a==b: ", a==b)
    print("a!=b: ", a!=b)
    print("a.is_zero(): ", a.is_zero())
    print("z.is_zero(): ", z.is_zero())
    last = None
    for k in blocks1.keys():
        last = k
        print("a("+str(k)+"): ", a.get_block(k))
    print("a("+str(last)+"): ", a(last))
    r = vec_class.zero_vec(blocks1.block_len(last))
    print("setting a("+str(last)+") to 0.")
    a.set_block(last,r)
    print("a = ",a)
    # testing with non-matching blocks
    print("=============================")
    a = BlockVec(blocks1, vec_class.rand(len(blocks1)))
    b = BlockVec(blocks2, vec_class.rand(len(blocks1)))
    print("a = ", a)
    print("b = ", b)
    print("b.clone() = ", b.clone())
    print("Setting a[1] = b[1] and a[2] = b[2]")
    a.set_entry(1,b.get_entry(1))
    a[2] = b[2]
    print("a = ", a)
    print("b = ", b)
    print("b[3:5] = ", b.get_segment(3,5))
    print("setting a[5:7] to b[3:5]")
    a.set_segment(5, b.get_segment(3,5))
    print("a =   ", a)
    print("b =   ", b)
    print("a+b = ", a+b)
    print("a-b = ", a-b)
    print("-a  = ", a)
    print("a*b = ", a*b)
    print("a==b: ", a==b)
    print("a!=b: ", a!=b)
    c = BlockVec(blocks2, a._vec)
    print("c: ", c)
    print("a==c: ", a==c)
    print("a!=c: ", a!=c)
    # testing with non-blocked vectors
    print("=============================")
    a = BlockVec(blocks1, vec_class.rand(len(blocks1)))
    b = vec_class.rand(len(blocks1))
    print("a = ", a)
    print("b = ", b)
    print("b[3:5] = ", b.get_segment(3,5))
    print("setting a[5:7] to b[3:5]")
    a.set_segment(5, b.get_segment(3,5))
    print("a =   ", a)
    print("b =   ", b)
    print("a+b = ", a+b)
    print("a-b = ", a-b)
    print("-a  = ", a)
    print("a*b = ", a*b)
    print("a==b: ", a==b)
    print("a!=b: ", a!=b)
    c = a._vec.clone()
    print("c: ", c)
    print("a==c: ", a==c)
    print("a!=c: ", a!=c)

def _test_BitVec(vec_class1, vec_class2):
    z = vec_class1.zero_vec(8)
    a = vec_class1.rand(8)
    b = vec_class2.rand(8)
    # Visual tests.
    # Arithmetic tests.
    print("a   = ", a)
    print("b   = ", b)
    print("a+b = ", a+b)
    print("a-b = ", a-b)
    print("-a  = ", -a)
    print("a*b = ", a*b)
    print("len(a)=", len(a))
    # Get/Set tests
    print("a[3] = ", a[3])
    print("b[7] = ", b[7])
    print("b[5] = ", b.get_entry(5))
    print("b[2:6] = ", b.get_segment(2,6))
    print("Setting b[2:6] to a[1:5].")
    b.set_segment(2, a.get_segment(1,5))
    print("a = ", a)
    print("b = ", b)
    print("Setting b[7] to a[6].")
    b[7] = a[6]
    print("a = ", a)
    print("b = ", b)
    print("Setting b[0] to a[7].")
    b.set_entry(0,a.get_entry(7))
    print("a = ", a)
    print("b = ", b)
    zz = a.vec_class().zero_vec(7)
    print("A zero vector of length 7:", zz)
    # Theoreatical tests.
    # Arithmetic tests
    res = []
    z = vec_class2.zero_vec(100)
    a = vec_class1.rand(100)
    b = vec_class1.rand(100)
    c = vec_class2.rand(100)
    res.append((1, z == a-a))
    res.append((2, not(z != b-b)))
    res.append((3, a+b == b+a))
    res.append((4, a*(b-c) == (a*b) ^ (a*c) ))
    res.append((5, not(a*z != z*b)))
    res.append((5.5, vec_class1.zero_vec(100) == z) )
    res.append((6, not(vec_class2.zero_vec(99) == z)))
    res.append((7, vec_class2.zero_vec(99) != z))
    res.append((8, z.is_zero()))
    res.append((9, not a.is_zero())) # might be False, but this is unlikely.
    res.append((10, not(a == b) )) # might be False, but this is unlikely.
    res.append((11, z == z.vec_class().zero_vec(100)))
    # get/set tests
    seg = a.get_segment(0,60)
    b.set_segment(0,seg)
    res.append((12, a.get_entry(25) == b.get_entry(25)))
    res.append((13, a.get_segment(13,19) == b.get_segment(13,19)))
    seg = b.get_segment(60,100)
    a.set_segment(60,seg)
    res.append((14, a == b))
    b.set_entry(14, b.get_entry(15))
    res.append((15, b[14] == b[15]))
    res.append((16, len(a)==100))
    res.append( (17, a.vec_class() == vec_class1) )
    return (False not in [y for x,y in res]), res

def _test_cohomology(mat_class):
    A = mat_class.zero_mat(20,29)
    B = mat_class.zero_mat(29,15)
    res = []
    for i in range(0,10):
        A.set_to_one(i,i)
    for j in range(0,10):
        B.set_to_one(j+10,j)
    X = mat_class.rand_inv_mat(20)
    Y = mat_class.rand_inv_mat(29)
    Yinv = Y.quasi_inverse()
    Z = mat_class.rand_inv_mat(15)
    #
    res.append( (1, Y*Yinv == mat_class.const_id_mat(29)) )
    res.append( (2, (A*B).is_zero() ) )
    d0 = Y*B*Z
    d1 = X*A*Yinv
    res.append( (2.5, (d1*d0).is_zero() ) )
    d0ker = d0.ker()
    d0coker = d0.coker()
    res.append( (3, d0ker.cols() == 5 and d0ker.rows() == 15) )
    res.append( (4, d0coker.rows() == 19 and d0coker.cols() == 29) )
    #
    coh1 = mat_class.cohomology_basis(d1,d0)
    coh2 = mat_class.rand_cohomology_basis(d1,d0)
    dual = mat_class.cohomology_dual_basis(d1,d0)
    #
    res.append( (5, (d1 * coh1.transpose()).is_zero() ) )
    res.append( (6, (coh1.rows() == 9) and (coh1.cols() == 29) ) )
    B = d0.transpose()
    B.append_rows(coh1)
    res.append( (7, B.Gauss_elim().rows() == 10 + coh1.rows() ) )
    #
    res.append( (8, (d1 * coh2.transpose()).is_zero() ) )
    res.append( (9, (coh2.rows() == 9) and (coh1.cols() == 29) ) )
    B = d0.transpose()
    B.append_rows(coh2)
    res.append( (10, B.Gauss_elim().rows() == 10 + coh2.rows() ) )
    res.append( (11, coh1 != coh2) )
    #
    res.append( (12, (dual*d0).is_zero() ) )
    res.append( (13, (dual.rows() == 9) and (dual.cols() == 29) ) )
    res.append( (14, (dual* coh1.transpose()).rank() == 9) )
    res.append( (15, (dual* coh2.transpose()).rank() == 9) )
    return (False not in [y for x,y in res]), [x for x,y in res if y == False] 

def _test_BitMat(mat_class1, mat_class2):
    Z = mat_class1.zero_mat(40,60)
    A = mat_class1.rand(40,60)
    B = mat_class2.rand(40,60)
    C = mat_class2.rand(60,40)
    res = []
    res.append( (1,Z.is_zero()) )
    res.append( (2, not A == B) )
    res.append( (3, A != B) )
    res.append( (4, A+B == B+A) )
    res.append( (5, -(B-A) == (A-B+Z)) )
    res.append( (6, (A+B)*C == (A*C) + (B*C)) )
    A.set_to_one(10,10)
    B.set_to_one(20,10)
    res.append( (7, A[10,10] == B[20, 10]) )
    A.set_to_zero(30, 30)
    B.set_to_zero(30,40)
    res.append( (8, A.get_entry(30,30) == B.get_entry(30, 40)) )
    A.set_entry(15,15,B.get_entry(14,14))
    res.append( (9, A[15,15] == B.get_entry(14,14)) )
    B_clone = B.clone()
    G1 = B.Gauss_elim()
    G2, l1 = B.Gauss_elim_extended()
    G3, l2, E = B.Gauss_elim_extended_extras(B)
    res.append( (10, B_clone == B) )
    res.append( (10.5, G1 == G2 and G2 == G3) )
    res.append( (11, l1 == l2) )
    res.append( (12, E == G3) )
    B.append_rows(B)
    G4, l3, Z2 = B.Gauss_elim_extended_extras(B.mat_class().zero_mat(40+40,60))
    res.append( (13, G1 == G4) )
    res.append( (14, l1 == l3) )
    res.append( (15, Z2.is_zero()) )
    B = B_clone
    B.add_zero_rows(30)
    G5, l4 = B.Gauss_elim_extended()
    res.append( (16, G1 == G5) )
    res.append( (17, l1 == l4) )
    C = B.mat_class().rand_inv_mat(61)
    I = B.mat_class().const_id_mat(61)
    res.append( (18, not I.is_zero()) )
    G, l, inv = C.Gauss_elim_extended_extras(I)
    res.append( (18, I == C.mat_class().id_mat(61)) )
    res.append( (19, inv * C == I) )
    res.append( (20, not G != I) )
    res.append( (21, len(l) == 61) )
    A = B.mat_class().rand(40,60)
    Q = A.quasi_inverse()
    AQ = A*Q
    res.append( (22, AQ == Q.mat_class().id_mat(40)) )
    res.append( (23, not AQ == Q.mat_class().const_id_mat(41)) )
    res.append( (24, Q*A != Q.mat_class().const_id_mat(60)) )
    res.append( (25, AQ*A == A) )
    C = C.mat_class().rand(60,40)
    Q = C.quasi_inverse()
    res.append( (26, C*Q*C == C) )
    D = A.mat_class().rand(30,30)
    E = D - D.transpose()
    res.append( (27, E == -E.transpose()) )
    res.append( (28, E.get_row(10) == E.get_col(10)) )
    N = A.null_space()
    Nt = N.transpose()
    res.append( (29, (A*Nt).is_zero()) )
    G = N.Gauss_elim()
    res.append( (30, G.rows() == N.rows() and N.rows() == Nt.cols() ) )
    A.append_rows(A) # now A is 80-by-60 of rank 40
    K = A.ker()
    C = A.coker()
    res.append( (31, (A*K).is_zero()) )
    res.append( (32, (C*A).is_zero()) )
    res.append( (33, K.cols() + C.rows() == A.cols() ) )
    v = K.get_col(1)
    res.append( (34, A(v).is_zero() ) )
    u = A.get_col(1)
    res.append( (35, C.apply(u).is_zero() ) )
    #
    A = mat_class1.zero_mat(21,21)
    A.set_to_one(0,0)
    A.set_to_one(11,13)
    res.append( (36, A[0,0] == A[11,13] ) )
    A.set_to_one(20,20)
    A.set_to_one(11,13)
    res.append( (37, A[0,0] == A[11,13] ) )
    A.set_to_zero(11,13)
    res.append( (38, A[19,20] == A[11,13] ) )
    A.set_to_zero(0,0)
    A.set_to_zero(5,0)
    res.append( (38, not A.is_zero() ) )
    #
    U = mat_class1.rand_inv_mat(20)
    A = mat_class2.zero_mat(40,50)
    A.set_rect(0,20,U)
    A.set_row_segment(20,20,U.get_row(19))
    res.append( (39, not A.is_zero() ) )
    res.append( (40, not A.get_row(19) != A.get_row(20) ) )
    A.set_row(19, A.get_row(25))
    res.append( (41, A.get_row(19) == A.get_row(30) ) )
    U = mat_class2.rand_inv_mat(20)
    A = mat_class2.zero_mat(50,50)
    A.set_rect(5,20,U)
    A.set_col_segment(5,0,U.get_col(1))
    A.set_col_segment(0,1,U.get_col(2))
    res.append( (42, A.get_col(0) == A.get_col(21) ) )
    res.append( (43, A.get_col(1) != A.get_col(22) ) )
    A.set_col(22, A.get_row(48))
    res.append( (44, A.get_col(22) == A.get_row(48) ) )
    res.append( (45, A.mat_class() == mat_class2) )
    res.append( (46, not A.mat_class() != mat_class2) )
    #
    A = mat_class1.rand(30,31)
    B = mat_class2.rand(15,16)
    A.set_rect(2,3,B)
    R = A.get_rect(5,11,5,12)
    res.append( (46.1, R == B.get_rect(5-2,11-2,5-3,12-3)) )
    res.append( (46.2, R.rows() == 11-5 and R.cols() == 12-5) )
    res.append( (46.3, R.mat_class() == mat_class1) )
    A = mat_class1.zero_mat(30,31)
    res.append( (46.4, A.get_rect(10,14,5,20).is_zero()) )
    #
    vec_class = mat_class1.vec_class()
    v = vec_class.rand(30)
    u = vec_class.rand(30)
    A = mat_class2.rand(25,30)
    res.append( (47, A(u+v) == A(u) + A.apply(v)) )
    res.append( (48, not A(v).is_zero() ) )
    #
    A = mat_class1.rand(25,30)
    res.append( (49, len(A.get_row(1)) == A.cols()) )
    res.append( (50, len(A.get_col(2)) == A.rows()) )
    res.append( (51, A.get_col(2).vec_class() == mat_class1.vec_class()) )
    # testing zero const matrix
    A = mat_class1.const_zero_mat(20,30)
    res.append( (51.5, A.is_zero() ) )
    res.append( (52, A.transpose() == mat_class2.const_zero_mat(30,20)) )
    res.append( (53, A.Gauss_elim().rows() == 0 ) )
    res.append( (53.5, A.Gauss_elim().cols() == A.cols() ) )
    res.append( (54, A.mat_class() == mat_class1) )
    B = mat_class2.rand(20,15)
    G, l, E = A.Gauss_elim_extended_extras(B)
    res.append( (55, B == E) )
    res.append( (56, G.rows() == 0 and G.cols() == A.cols() ) )
    res.append( (57, len(l) == 0) )
    G, l = A.Gauss_elim_extended()
    res.append( (58, G.rows() == 0 and G.cols() == A.cols() ) )
    res.append( (59, len(l) == 0) )
    res.append( (60, A.inverse() is None) )
    B = mat_class2.rand(20,30)
    res.append( (61, A+B == B and B+A == B) )
    res.append( (62, A-B == -B and B-A == B) )
    res.append( (63, A * mat_class1.rand(30,11) == mat_class2.zero_mat(20,11) ) )
    res.append( (64, -A == A) )
    res.append( (65, mat_class1.rand(11,20) * A == mat_class2.const_zero_mat(11,30) ) )
    v = mat_class2.vec_class().rand(30)
    res.append( (65.1, A(v).is_zero() ) )
    # testing id const matrix
    A = mat_class1.const_id_mat(22)
    res.append( (65.5, not A.is_zero() ) )
    res.append( (66, A.transpose() == mat_class2.const_id_mat(22)) )
    res.append( (67, A.Gauss_elim().rows() == 22 ) )
    res.append( (68, A.Gauss_elim().cols() == A.cols() ) )
    res.append( (69, A.mat_class() == mat_class1) )
    B = mat_class2.rand(22,15)
    G, l, E = A.Gauss_elim_extended_extras(B)
    res.append( (70, B == E) )
    res.append( (71, G.rows() == A.rows() and G.cols() == A.cols() ) )
    res.append( (72, len(l) == A.rows()) )
    G, l = A.Gauss_elim_extended()
    res.append( (73, G.rows() == 22 and G.cols() == 22 ) )
    res.append( (74, len(l) == 22) )
    res.append( (75, A.inverse() == mat_class1.id_mat(22)) )
    B = mat_class2.rand(22,22)
    res.append( (76, A+B == B+A) )
    res.append( (77, A-B == -(B-A)) )
    B = mat_class2.rand(22,31)
    C = mat_class2.rand(23,22)
    res.append( (78, A * B == B ) )
    res.append( (79, C * A == C) )
    res.append( (80, mat_class1.const_zero_mat(17,22) * A == mat_class2.zero_mat(17,22) ) )
    v = mat_class2.vec_class().rand(22)
    res.append( (80.1, A(v) == v ) )
    # test __eq__ and __neq__ for zero and identity matrix
    res.append( (80.11, mat_class1.const_zero_mat(5,6) == mat_class2.const_zero_mat(5,6)) )
    res.append( (80.12, not mat_class1.const_zero_mat(5,6) == mat_class2.const_zero_mat(5,5)) )
    res.append( (80.13, not mat_class1.const_zero_mat(5,6) == mat_class2.const_id_mat(5)) )
    res.append( (80.14, not mat_class1.const_zero_mat(5,6) == mat_class2.id_mat(5)) )
    res.append( (80.15, not mat_class1.const_zero_mat(5,6) == mat_class2.rand(5,6)) )
    #
    res.append( (80.16, not mat_class1.const_zero_mat(5,6) != mat_class2.const_zero_mat(5,6)) )
    res.append( (80.17, mat_class1.const_zero_mat(5,6) != mat_class2.const_zero_mat(5,5)) )
    res.append( (80.18, mat_class1.const_zero_mat(5,6) != mat_class2.const_id_mat(5)) )
    res.append( (80.19, mat_class1.const_zero_mat(5,6) != mat_class2.id_mat(5)) )
    res.append( (80.20, mat_class1.const_zero_mat(5,6) != mat_class2.rand(5,6)) )
    #
    res.append( (80.21, mat_class1.const_id_mat(5) == mat_class2.const_id_mat(5)) )
    res.append( (80.22, not mat_class1.const_id_mat(5) == mat_class2.id_mat(6)) )
    res.append( (80.23, not mat_class1.const_id_mat(5) == mat_class2.zero_mat(5,5)) )
    res.append( (80.24, not mat_class1.const_id_mat(5) == mat_class2.rand(5,5)) )
    res.append( (80.25, not mat_class1.const_id_mat(5) == mat_class2.const_zero_mat(5,5)) )
    #
    res.append( (80.26, not mat_class1.const_id_mat(5) != mat_class2.const_id_mat(5)) )
    res.append( (80.27, mat_class1.const_id_mat(5) != mat_class2.id_mat(6)) )
    res.append( (80.28, mat_class1.const_id_mat(5) != mat_class2.zero_mat(5,5)) )
    res.append( (80.29, mat_class1.const_id_mat(5) != mat_class2.rand(5,5)) )
    res.append( (80.30, mat_class1.const_id_mat(5) != mat_class2.const_zero_mat(5,5)) )
    #
    res.append( (80.31, mat_class1.const_id_mat(0) == mat_class2.const_zero_mat(0,0)) )
    res.append( (80.32, not mat_class1.const_id_mat(0) != mat_class2.const_zero_mat(0,0)) )
    res.append( (80.33, mat_class1.const_zero_mat(0,0) == mat_class2.const_id_mat(0)) )
    res.append( (80.34, not mat_class1.const_zero_mat(0,0) != mat_class2.const_id_mat(0)) )
    # test mutablitly
    A = mat_class1.const_zero_mat(11,12).to_mutable()
    A.set_to_one(10,11) # all  is good if it didn't collapse
    res.append( (81, A.is_mutable() ) )
    res.append( (82, A.to_mutable() is A) )
    res.append( (83, A.mutable_clone() is not A) )
    A = mat_class1.const_id_mat(12).to_mutable()
    A.set_row_segment(11,0,mat_class2.vec_class().zero_vec(12)) # all  is good if it didn't collapse
    res.append( (84, A.is_mutable() ) )
    res.append( (85, A.to_mutable() is A) )
    res.append( (86, A.mutable_clone() is not A) )
    # test get_row_as, get_col_as
    M = mat_class1.rand(20,22)
    res.append( (87, M.get_row(5) == M.get_row_as(5,mat_class2.vec_class())) )
    res.append( (88, M.get_col(5) == M.get_col_as(5,mat_class2.vec_class())) )
    res.append( (89, M.get_row_as(5,mat_class2.vec_class()).vec_class() == mat_class2.vec_class()) )
    res.append( (90, M.get_col_as(5,mat_class2.vec_class()).vec_class() == mat_class2.vec_class()) )
    return (False not in [y for x,y in res]), [x for x,y in res if y == False] 


def _test_BaseVector(vec_class):
    z = vec_class.zero_vec(8)
    a = vec_class.rand(8)
    b = vec_class.rand(8)
    # Visual tests.
    # Arithmetic tests.
    print("a   = ", a)
    print("b   = ", b)
    print("a+b = ", a+b)
    print("a-b = ", a-b)
    print("-a  = ", -a)
    print("a*b = ", a*b)
    print("len(a)=", len(a))
    # Get/Set tests
    print("a[3] = ", a[3])
    print("b[7] = ", b[7])
    print("b[5] = ", b.get_entry(5))
    print("b[2:6] = ", b.get_segment(2,6))
    print("Setting b[2:6] to a[1:5].")
    b.set_segment(2, a.get_segment(1,5))
    print("a = ", a)
    print("b = ", b)
    print("Setting b[7] to a[6].")
    b[7] = a[6]
    print("a = ", a)
    print("b = ", b)
    print("Setting b[0] to a[7].")
    b.set_entry(0,a.get_entry(7))
    print("a = ", a)
    print("b = ", b)
    zz = a.vec_class().zero_vec(7)
    print("A zero vector of length 7:", zz)
    # Theoreatical tests.
    # Arithmetic tests
    res = []
    z = vec_class.zero_vec(100)
    a = vec_class.rand(100)
    b = vec_class.rand(100)
    c = vec_class.rand(100)
    res.append((1, z == a-a))
    res.append((2, not(z != b-b)))
    res.append((3, a+b == b+a))
    res.append((4, a*(b-c) == a*b - a*c))
    res.append((5, not(a*z != z*b)))
    res.append((6, not(vec_class.zero_vec(99) == z)))
    res.append((7, vec_class.zero_vec(99) != z))
    res.append((8, z.is_zero()))
    res.append((9, not a.is_zero())) # might be False, but this is unlikely.
    res.append((10, not(a == b) )) # might be False, but this is unlikely.
    res.append((11, z == z.vec_class().zero_vec(100)))
    # get/set tests
    seg = a.get_segment(0,60)
    b.set_segment(0,seg)
    res.append((12, a.get_entry(25) == b.get_entry(25)))
    res.append((13, a.get_segment(13,19) == b.get_segment(13,19)))
    seg = b.get_segment(60,100)
    a.set_segment(60,seg)
    res.append((14, a == b))
    b.set_entry(14, b.get_entry(15))
    res.append((15, b[14] == b[15]))
    res.append((16, len(a)==100))
    res.append((17, a.vec_class() == vec_class) )
    res.append((18, not a.vec_class() != vec_class) )
    v = vec_class.rand(20)
    res.append((19, v == v.clone()) )
    # test scale
    a = vec_class.rand(100)
    a_clone = a.clone()
    b = vec_class.rand(100)
    s = vec_class.rand(1)[0] # a fishy way to get a random scalar.
    t = vec_class.rand(1)[0]
    res.append( (20, a.scale(s) + b.scale(s) == (a+b).scale(s)) )
    res.append( (21, a.scale(t) - b.scale(t) == (a-b).scale(t)) )
    res.append( (22, a.scale(s).scale(t) == a.scale(s*t)) )
    res.append( (23, a == a_clone) )
    return (False not in [y for x,y in res]), [x for x,y in res if y == False]




def _test_BaseMatrix(mat_class):
    Z = mat_class.zero_mat(40,60)
    A = mat_class.rand(40,60)
    B = mat_class.rand(40,60)
    C = mat_class.rand(60,40)
    res = []
    #
    res.append( (1,Z.is_zero()) )
    res.append( (2, not A == B) )
    res.append( (3, A != B) )
    res.append( (4, A+B == B+A) )
    res.append( (5, -(B-A) == (A-B+Z)) )
    res.append( (6, (A+B)*C == (A*C) + (B*C)) )
    A.set_to_one(10,10)
    B.set_to_one(20,10)
    res.append( (7, A[10,10] == B[20, 10]) )
    A.set_to_zero(30, 30)
    B.set_to_zero(30,40)
    res.append( (8, A.get_entry(30,30) == B.get_entry(30, 40)) )
    A.set_entry(15,15,B.get_entry(14,14))
    res.append( (9, A[15,15] == B.get_entry(14,14)) )
    #
    B_clone = B.clone()
    G1 = B.Gauss_elim()
    G2, l1 = B.Gauss_elim_extended()
    G3, l2, E = B.Gauss_elim_extended_extras(B)
    res.append( (10, B_clone == B) )
    res.append( (10.5, G1 == G2 and G2 == G3) )
    res.append( (11, l1 == l2) )
    res.append( (12, E == G3) )
    B.append_rows(B)
    G4, l3, Z2 = B.Gauss_elim_extended_extras(B.mat_class().zero_mat(40+40,60))
    res.append( (13, G1 == G4) )
    res.append( (14, l1 == l3) )
    res.append( (15, Z2.is_zero()) )
    B = B_clone
    B.add_zero_rows(30)
    G5, l4 = B.Gauss_elim_extended()
    res.append( (16, G1 == G5) )
    res.append( (17, l1 == l4) )
    C = B.mat_class().rand_inv_mat(61)
    I = B.mat_class().const_id_mat(61)
    res.append( (18, not I.is_zero()) )
    G, l, inv = C.Gauss_elim_extended_extras(I)
    res.append( (18, I == C.mat_class().id_mat(61)) )
    res.append( (19, inv * C == I) )
    res.append( (20, not G != I) )
    res.append( (21, len(l) == 61) )
    A = B.mat_class().rand(40,60)
    Q = A.quasi_inverse()
    AQ = A*Q
    res.append( (22, AQ == Q.mat_class().id_mat(40)) )
    res.append( (23, not AQ == Q.mat_class().const_id_mat(41)) )
    res.append( (24, Q*A != Q.mat_class().const_id_mat(60)) )
    res.append( (25, AQ*A == A) )
    C = C.mat_class().rand(60,40)
    Q = C.quasi_inverse()
    res.append( (26, C*Q*C == C) )
    D = A.mat_class().rand(30,30)
    E = D - D.transpose()
    res.append( (27, E == -E.transpose()) )
    res.append( (28, E.get_row(10) == -E.get_col(10)) )
    N = A.null_space()
    Nt = N.transpose()
    res.append( (29, (A*Nt).is_zero()) )
    G = N.Gauss_elim()
    res.append( (30, G.rows() == N.rows() and N.rows() == Nt.cols() ) )
    #
    A.append_rows(A) # now A is 80-by-60 of rank 40
    K = A.ker()
    C = A.coker()
    res.append( (31, (A*K).is_zero()) )
    res.append( (32, (C*A).is_zero()) )
    res.append( (33, K.cols() + C.rows() == A.cols() ) )
    v = K.get_col(1)
    res.append( (34, A(v).is_zero() ) )
    u = A.get_col(1)
    res.append( (35, C.apply(u).is_zero() ) )
    #
    A = mat_class.zero_mat(21,21)
    A.set_to_one(0,0)
    A.set_to_one(11,13)
    res.append( (36, A[0,0] == A[11,13] ) )
    A.set_to_one(20,20)
    A.set_to_one(11,13)
    res.append( (37, A[0,0] == A[11,13] ) )
    A.set_to_zero(11,13)
    res.append( (38, A[19,20] == A[11,13] ) )
    A.set_to_zero(0,0)
    A.set_to_zero(5,0)
    res.append( (38, not A.is_zero() ) )
    # 
    U = A.mat_class().rand_inv_mat(20)
    A = mat_class.zero_mat(40,50)
    A.set_rect(0,20,U)
    A.set_row_segment(20,20,U.get_row(19))
    res.append( (39, not A.is_zero() ) )
    res.append( (40, not A.get_row(19) != A.get_row(20) ) )
    A.set_row(19, A.get_row(25))
    res.append( (41, A.get_row(19) == A.get_row(30) ) )
    U = A.mat_class().rand_inv_mat(20)
    A = mat_class.zero_mat(50,50)
    A.set_rect(5,20,U)
    A.set_col_segment(5,0,U.get_col(1))
    A.set_col_segment(0,1,U.get_col(2))
    res.append( (42, A.get_col(0) == A.get_col(21) ) )
    res.append( (43, A.get_col(1) != A.get_col(22) ) )
    A.set_col(22, A.get_row(48))
    res.append( (44, A.get_col(22) == A.get_row(48) ) )
    res.append( (45, A.mat_class() == mat_class) )
    res.append( (46, not A.mat_class() != mat_class) )
    #
    A = mat_class.rand(30,31)
    B = mat_class.rand(15,16)
    A.set_rect(2,3,B)
    R = A.get_rect(5,11,5,12)
    res.append( (46.1, R == B.get_rect(5-2,11-2,5-3,12-3)) )
    res.append( (46.2, R.rows() == 11-5 and R.cols() == 12-5) )
    res.append( (46.3, R.mat_class() == mat_class) )
    A = mat_class.zero_mat(30,31)
    res.append( (46.4, A.get_rect(10,14,5,20).is_zero()) )
    #
    vec_class = mat_class.vec_class()
    v = vec_class.rand(30)
    u = vec_class.rand(30)
    A = mat_class.rand(25,30)
    res.append( (47, A(u+v) == A(u) + A.apply(v)) )
    res.append( (48, not A(v).is_zero() ) )
    # test get_row and get_col's length and vector class
    A = mat_class.rand(25,30)
    res.append( (49, len(A.get_row(1)) == A.cols()) )
    res.append( (50, len(A.get_col(2)) == A.rows()) )
    res.append( (51, A.get_col(2).vec_class() == mat_class.vec_class()) )
    # testing zero const matrix
    A = mat_class.const_zero_mat(20,30)
    res.append( (52, A.transpose() == mat_class.const_zero_mat(30,20)) )
    res.append( (53, A.Gauss_elim().rows() == 0 ) )
    res.append( (53.5, A.Gauss_elim().cols() == A.cols() ) )
    res.append( (54, A.mat_class() == mat_class) )
    B = mat_class.rand(20,15)
    G, l, E = A.Gauss_elim_extended_extras(B)
    res.append( (55, B == E) )
    res.append( (56, G.rows() == 0 and G.cols() == A.cols() ) )
    res.append( (57, len(l) == 0) )
    G, l = A.Gauss_elim_extended()
    res.append( (58, G.rows() == 0 and G.cols() == A.cols() ) )
    res.append( (59, len(l) == 0) )
    res.append( (60, A.inverse() is None) )
    B = mat_class.rand(20,30)
    res.append( (61, A+B == B and B+A == B) )
    res.append( (62, A-B == -B and B-A == B) )
    res.append( (63, A * mat_class.rand(30,11) == mat_class.zero_mat(20,11) ) )
    res.append( (64, -A == A) )
    res.append( (65, mat_class.rand(11,20) * A == mat_class.const_zero_mat(11,30) ) )
    # testing id const matrix
    A = mat_class.const_id_mat(22)
    res.append( (66, A.transpose() == mat_class.const_id_mat(22)) )
    res.append( (67, A.Gauss_elim().rows() == 22 ) )
    res.append( (68, A.Gauss_elim().cols() == A.cols() ) )
    res.append( (69, A.mat_class() == mat_class) )
    B = mat_class.rand(22,15)
    G, l, E = A.Gauss_elim_extended_extras(B)
    res.append( (70, B == E) )
    res.append( (71, G.rows() == A.rows() and G.cols() == A.cols() ) )
    res.append( (72, len(l) == A.rows()) )
    G, l = A.Gauss_elim_extended()
    res.append( (73, G.rows() == 22 and G.cols() == 22 ) )
    res.append( (74, len(l) == 22) )
    res.append( (75, A.inverse() == mat_class.id_mat(22)) )
    B = mat_class.rand(22,22)
    res.append( (76, A+B == B+A) )
    res.append( (77, A-B == -(B-A)) )
    B = mat_class.rand(22,31)
    C = mat_class.rand(23,22)
    res.append( (78, A * B == B ) )
    res.append( (79, C * A == C) )
    res.append( (80, mat_class.const_zero_mat(17,22) * A == mat_class.zero_mat(17,22) ) )
    # test mutablitly
    A = mat_class.const_zero_mat(11,12).to_mutable()
    A.set_to_one(10,11) # all  is good if it didn't collapse
    res.append( (81, A.is_mutable() ) )
    res.append( (82, A.to_mutable() is A) )
    res.append( (83, A.mutable_clone() is not A) )
    A = mat_class.const_id_mat(12).to_mutable()
    A.set_row_segment(11,0,mat_class.vec_class().zero_vec(12)) # all  is good if it didn't collapse
    res.append( (84, A.is_mutable() ) )
    res.append( (85, A.to_mutable() is A) )
    res.append( (86, A.mutable_clone() is not A) )
    return (False not in [y for x,y in res]), [x for x,y in res if not y] 

def _test_BaseMatrix_and_BaseVector_visually(mat_class, vec_class):
    A = mat_class.rand_mat(3,4)
    B = mat_class.rand_mat(4,5)
    C = mat_class.rand_mat(3,4)
    print("A = ", A)
    print("B = ", B)
    print("C = ", C)
    print("A + C = ", A+C)
    print("A - C = ", A-C)
    print("A * B = ", A*B)
    print("-A = ", -A)
    print("Gauss elimination of A:", A.Gauss_elim())
    print("Gauss elimination of B:", B.Gauss_elim())
    print("Transpose of A:", A.transpose())
    print("A == B:", A==B)
    print("A == C:", A==C)
    print("A != B:", A!=B)
    print("A != C:", A!=C)
    print("A is zero:", A.is_zero())
    print("null space of A:", A.null_space())
    print("ker(A):", A.ker())
    print("coker(A):", A.coker())

def _test_visually_ConstBaseMatrix(mat_class):
    Z = mat_class.zero_mat(4,6)
    A = mat_class.rand(4,6)
    B = mat_class.rand(4,6)
    C = mat_class.rand(6,4)
    print("Z = ", Z)
    print("A = ", A)
    print("B = ", B)
    print("C = ", C)
    print("A+B = ", A+B)
    print("A-B = ", A-B)
    print("A*C = ", A*C)
    print("-A = ", -A)
    print("Gauss elimination of B:", B.Gauss_elim())
    G, leading = B.Gauss_elim_extended()
    print("Exteded Gauss elimination of B:", G,   leading)
    print("B = ", B)
    G, leading, M = B.Gauss_elim_extended_extras(A)
    print("Perform Gauss elimination on B and same operations on A:")
    print(G)
    print(M)

def _test_Matrix(p=5):
    F = PrimeField(p)
    return _test_BaseMatrix(MatrixClass(F))

def _test_Vector(p=5):
    F = PrimeField(p)
    return _test_BaseVector(VectorClass(F))
