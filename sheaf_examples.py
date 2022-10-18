# This file includes examples of sheaves that were investigated in a paper by First and Kaufman
# each function returns a sheaf and is named as
# <base_complex>_sheafdim<sheaf_dim>_<dims of 0th, 1st, 2nd cohomology>_<nickname>.

from sheaf_tools import StructureSheaf, DirectSumSheaf
from bit_mat import BitMat, BitMatClass
from simplicial_complex import torus_ndim_morphism
from ramanujan_complex import LSV_quo_morphism_2_to_1


def torus2dim_sheafdim41_1_2_1():
    cov = torus_ndim_morphism(grid=(3,3), factors=(3,3))
    X = cov.target
    F = StructureSheaf(X, BitMatClass(240))
    unit = F.unit_map(cov)
    coker = unit.coker()
    G = coker.target
    return DirectSumSheaf([F, G, G, G, G, G])

def torus3dim_sheafdim27_1_3_3():
    cov = torus_ndim_morphism(grid=(3,3,3), factors=(3,3,3))
    Y = cov.source
    F = StructureSheaf(Y, BitMatClass(240))
    return F.pushforward(cov)

def torus3dim_sheafdim40_0_0_0():
    cov = torus_ndim_morphism(grid=(3,3,3), factors=(1,3,3))
    X = cov.target
    F = StructureSheaf(X, BitMatClass(240))
    counit = F.counit_map(cov)
    ker = counit.ker()
    G = ker.source # dim G = 8
    return DirectSumSheaf([G, G, G, G, G])

def torus3dim_sheafdim49_1_3_3_mixA():
    cov1 = torus_ndim_morphism(grid=(3,3,3), factors=(1,3,3))
    cov2 = torus_ndim_morphism(grid=(3,3,3), factors=(3,1,3))
    cov3 = torus_ndim_morphism(grid=(3,3,3), factors=(3,3,1))
    X = cov1.target
    cov2.target = X
    cov3.target = X
    F = StructureSheaf(X, BitMatClass(240))
    counit1 = F.counit_map(cov1)
    ker1 = counit1.ker()
    G1 = ker1.source # dim G1 = 8
    counit2 = F.counit_map(cov2)
    ker2 = counit2.ker()
    G2 = ker2.source # dim G2 = 8
    counit3 = F.counit_map(cov3)
    ker3 = counit3.ker()
    G3 = ker3.source # dim G3 = 8
    return DirectSumSheaf([F, G1, G1, G2, G2, G3, G3])

def torus3dim_sheafdim37_1_3_3_mixB():
    cov1 = torus_ndim_morphism(grid=(3,3,3), factors=(1,3,3))
    cov2 = torus_ndim_morphism(grid=(3,3,3), factors=(5,1,1))
    cov3 = torus_ndim_morphism(grid=(3,3,3), factors=(1,1,7))
    X = cov1.target
    cov2.target = X
    cov3.target = X
    F = StructureSheaf(X, BitMatClass(240))
    counit1 = F.counit_map(cov1)
    ker1 = counit1.ker()
    G1 = ker1.source # dim G1 = 8
    counit2 = F.counit_map(cov2)
    ker2 = counit2.ker()
    G2 = ker2.source # dim G2 = 4
    counit3 = F.counit_map(cov3)
    ker3 = counit3.ker()
    G3 = ker3.source # dim G3 = 6
    return DirectSumSheaf([F, G1, G1, G2, G2, G3, G3])

def torus3dim_sheafdim42_2_6_6():
    cov = torus_ndim_morphism(grid=(3,3,3), factors=(3,3,1))
    X = cov.target
    F = StructureSheaf(X, BitMatClass(240))
    unit = F.unit_map(cov)
    coker = unit.coker()
    G = coker.target # dim(G) = 8
    return DirectSumSheaf([F, F, G, G, G, G, G])

def torus3dim_sheafdim35_1_4_5():
    cov1 = torus_ndim_morphism(grid=(3,3,3), factors=(3,3,1))
    cov2 = torus_ndim_morphism(grid=(3,3,3), factors=(2,2,3))
    cov2.target = cov1.target
    X = cov1.target
    F = StructureSheaf(X, BitMatClass(240))
    unit = F.unit_map(cov1)
    coker = unit.coker()
    G1 = coker.target # dim(G1) = 8
    counit = F.counit_map(cov2)
    ker = counit.ker()
    G2 = ker.source # dim(G2) = 11
    return DirectSumSheaf([G2, G1, G1, G1])

def torus3dim_sheafdim46_2_8_10():
    cov1 = torus_ndim_morphism(grid=(3,3,3), factors=(3,3,1))
    cov2 = torus_ndim_morphism(grid=(3,3,3), factors=(2,2,3))
    cov2.target = cov1.target
    X = cov1.target
    F = StructureSheaf(X, BitMatClass(240))
    unit = F.unit_map(cov1)
    coker = unit.coker()
    G1 = coker.target # dim(G1) = 8
    counit = F.counit_map(cov2)
    ker = counit.ker()
    G2 = ker.source # dim(G2) = 11
    return DirectSumSheaf([G2, G2, G1, G1, G1])

def torus3dim_sheafdim46_2_8_10b():
    cov1 = torus_ndim_morphism(grid=(3,3,3), factors=(3,1,3))
    cov2 = torus_ndim_morphism(grid=(3,3,3), factors=(2,2,3))
    cov2.target = cov1.target
    X = cov1.target
    F = StructureSheaf(X, BitMatClass(240))
    unit = F.unit_map(cov1)
    coker = unit.coker()
    G1 = coker.target # dim(G1) = 8
    counit = F.counit_map(cov2)
    ker = counit.ker()
    G2 = ker.source # dim(G2) = 11
    return DirectSumSheaf([G2, G2, G1, G1, G1])

def torus3dim_sheafdim34_1_3_3_even_partial():
    cov1 = torus_ndim_morphism(grid=(3,3,3), factors=(3,3,1))
    cov2 = torus_ndim_morphism(grid=(3,3,3), factors=(1,1,2))
    cov2.target = cov1.target
    X = cov1.target
    F1 = StructureSheaf(X, BitMatClass(240))
    unit = F1.unit_map(cov1)
    coker = unit.coker()
    G1 = coker.target # dim(G1) = 8
    F2 = StructureSheaf(cov2.source, BitMatClass(240))
    G2 = F2.pushforward(cov2) # dim(G2) = 2
    return DirectSumSheaf([G2, G1, G1, G1, G1])

def torus3dim_sheafdim34_1_3_3_even_partial_2skel():
    cov1 = torus_ndim_morphism(grid=(3,3,3), factors=(3,3,1), skeleton=2)
    cov2 = torus_ndim_morphism(grid=(3,3,3), factors=(1,1,2), skeleton=2)
    cov2.target = cov1.target
    X = cov1.target
    F1 = StructureSheaf(X, BitMatClass(240))
    unit = F1.unit_map(cov1)
    coker = unit.coker()
    G1 = coker.target # dim(G1) = 8
    F2 = StructureSheaf(cov2.source, BitMatClass(240))
    G2 = F2.pushforward(cov2) # dim(G2) = 2
    return DirectSumSheaf([G2, G1, G1, G1, G1])

def torus3dim_sheafdim36_1_3_3_even_partial():
    cov1 = torus_ndim_morphism(grid=(3,3,3), factors=(3,3,1))
    cov2 = torus_ndim_morphism(grid=(3,3,3), factors=(1,2,2))
    cov2.target = cov1.target
    X = cov1.target
    F1 = StructureSheaf(X, BitMatClass(240))
    unit = F1.unit_map(cov1)
    coker = unit.coker()
    G1 = coker.target # dim(G1) = 8
    F2 = StructureSheaf(cov2.source, BitMatClass(240))
    G2 = F2.pushforward(cov2) # dim(G2) = 4
    return DirectSumSheaf([G2, G1, G1, G1, G1])

def torus3dim_sheafdim40_1_3_3_even():
    cov1 = torus_ndim_morphism(grid=(3,3,3), factors=(3,3,1))
    cov2 = torus_ndim_morphism(grid=(3,3,3), factors=(2,2,2))
    cov2.target = cov1.target
    X = cov1.target
    F1 = StructureSheaf(X, BitMatClass(240))
    unit = F1.unit_map(cov1)
    coker = unit.coker()
    G1 = coker.target # dim(G1) = 8
    F2 = StructureSheaf(cov2.source, BitMatClass(240))
    G2 = F2.pushforward(cov2) # dim(G2) = 8
    return DirectSumSheaf([G2, G1, G1, G1, G1])

def torus3dim_sheafdim48_1_3_3_even():
    cov1 = torus_ndim_morphism(grid=(3,3,3), factors=(3,3,1))
    cov2 = torus_ndim_morphism(grid=(3,3,3), factors=(2,2,2))
    cov2.target = cov1.target
    X = cov1.target
    F1 = StructureSheaf(X, BitMatClass(240))
    unit = F1.unit_map(cov1)
    coker = unit.coker()
    G1 = coker.target # dim(G1) = 8
    F2 = StructureSheaf(cov2.source, BitMatClass(240))
    G2 = F2.pushforward(cov2) # dim(G2) = 8
    return DirectSumSheaf([G2, G1, G1, G1, G1, G1])

def LSVq1_sheafdim17_1_3_4643():
    cov = LSV_quo_morphism_2_to_1()
    X = cov.target
    Y = cov.source
    F = StructureSheaf(X, BitMatClass(480))
    unit = F.unit_map(cov)
    coker = unit.coker()
    G = coker.target
    return DirectSumSheaf([F, G])

def LSVq1_sheafdim33_1_3_9011():
    cov = LSV_quo_morphism_2_to_1()
    X = cov.target
    Y = cov.source
    F = StructureSheaf(X, BitMatClass(480))
    unit = F.unit_map(cov)
    coker = unit.coker()
    G = coker.target
    return DirectSumSheaf([F, G, G])

def LSVq1_sheafdim49_1_3_13379():
    cov = LSV_quo_morphism_2_to_1()
    X = cov.target
    Y = cov.source
    F = StructureSheaf(X, BitMatClass(480))
    unit = F.unit_map(cov)
    coker = unit.coker()
    G = coker.target
    return DirectSumSheaf([F, G, G, G])

def LSVq1_sheafdim65_1_3_17747():
    cov = LSV_quo_morphism_2_to_1()
    X = cov.target
    Y = cov.source
    F = StructureSheaf(X, BitMatClass(480))
    unit = F.unit_map(cov)
    coker = unit.coker()
    G = coker.target
    return DirectSumSheaf([F, G, G, G, G])

def LSVq1_sheafdim50_2_6_9286():
    cov = LSV_quo_morphism_2_to_1()
    X = cov.target
    Y = cov.source
    F = StructureSheaf(X, BitMatClass(480))
    unit = F.unit_map(cov)
    coker = unit.coker()
    G = coker.target
    return DirectSumSheaf([F, F, G, G, G])

def LSVq1_sheafdim32_0_0_8736():
    cov = LSV_quo_morphism_2_to_1()
    X = cov.target
    Y = cov.source
    F = StructureSheaf(X, BitMatClass(480))
    unit = F.unit_map(cov)
    coker = unit.coker()
    G = coker.target
    return DirectSumSheaf([G, G])

def LSVq1_sheafdim80_0_0_21840():
    cov = LSV_quo_morphism_2_to_1()
    X = cov.target
    Y = cov.source
    F = StructureSheaf(X, BitMatClass(480))
    unit = F.unit_map(cov)
    coker = unit.coker()
    G = coker.target
    return DirectSumSheaf([G, G, G, G, G])
