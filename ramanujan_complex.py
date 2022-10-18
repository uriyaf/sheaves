# implements several quotients of the explicit Ramanujan complex in the paper by Lubotsky, Samuels and Vishne (LSV, 2005). 
# The main functions are:
# (*) LSV_quotient1, LSV_quotient2, LSVquotient3.
#   They return Ramanujan complexes with 273, 273*17, 273*17*3 vertices,
#   and each complex covers the one prior to it.
# (*) LSV_quo_morphism_2_to_1, LSV_quo_morphism_3_to_1.
#   They return the morphism from LSV_quotient<2 or 3> to LSV_quotient1.

from matrix import Vector
from simplicial_complex import DeltaSet, DeltaSetMorphism, Schrier_complex
from fields import F16
from matrix import Matrix, Vector, MatrixClass
from matrix import *


def proj_plane(field):
    """
    Returns iterable enumerating all points in the projective plane over the finite field 'field'.
    """
    one = field.one()
    zero = field.zero()
    for a in field.enum():
        for b in field.enum():
            yield (one,a,b)
    for a in field.enum():
        yield (zero,one,a)
    yield (zero,zero,one)

def normalize_vec(vec):
    """
    Takes a vector (an instance of Vector) and returns a scaling of the vector
    so that the first (left-most) nonzero entry is 1.
    """
    field = vec.get_field()
    zero = field.zero()
    i = 0
    while i < len(vec):
        if vec[i] != zero:
            break
        i += 1
    if i == len(vec):
        return vec.clone()
    a = vec[i].inverse()
    return Vector(field, vec=(a*vi for vi in vec))

def mul_mat_by_proj_vec(mat, vec):
    return normalize_vec(mat.apply(vec))

def is_mat_has_no_stable_proj_vecs(mat):
    """
    Takes a 3x3 matrix M over field. Returns True if there is no nonzero v such that Mv is proportional to v.
    """
    assert mat.rows() == 3 and mat.cols() == 3
    for vec in proj_plane(field):
        if vec == mul_mat_by_proj_vec(mat,vec):
            return False
    return True

def LSV_lattice_generators():
    b0 = Matrix(F16, 3, 3, mat=[[F16(2+8), F16(4), F16(2+4)], \
      [F16(2), F16(8), F16(1+2+4)], \
      [F16(2+4), F16(1+4), F16(1+8)]])
    b1 = Matrix(F16, 3, 3, mat=[[F16(1+2+4+8), F16(2+4), F16(1+4)], \
      [F16(1+2), F16(4+8), F16(1)], \
      [F16(1+4), F16(2), F16(8)]])
    b2 = Matrix(F16, 3, 3, mat=[[F16(1+4+8), F16(1+4), F16(2)], \
      [F16(1+2+4), F16(2+8), F16(4)], \
      [F16(2), F16(1+2), F16(4+8)]])
    b3 = Matrix(F16, 3, 3, mat=[[F16(2+4+8), F16(2), F16(1+2)], \
      [F16(1), F16(1+2+4+8), F16(2+4)], \
      [F16(1+2), F16(1+2+4), F16(2+8)]])
    b4 = Matrix(F16, 3, 3, mat=[[F16(1+8), F16(1+2), F16(1+2+4)], \
      [F16(4), F16(1+4+8), F16(1+4)], \
      [F16(1+2+4), F16(1), F16(1+2+4+8)]])
    b5 = Matrix(F16, 3, 3, mat=[[F16(8), F16(1+2+4), F16(1)], \
      [F16(2+4), F16(2+4+8), F16(2)], \
      [F16(1), F16(4), F16(1+4+8)]])
    b6 = Matrix(F16, 3, 3, mat=[[F16(4+8), F16(1), F16(4)], \
      [F16(1+4), F16(1+8), F16(1+2)], \
      [F16(4), F16(2+4), F16(2+4+8)]])
    return [b0,b1,b2,b3,b4,b5,b6]

def LSV_face_generators_and_labels():
    id3 = MatrixClass(F16).const_id_mat(3)  
    b = LSV_lattice_generators()
    gens = [[(id3,)], \
            [(id3, bi) for bi in b], \
            [[id3, b[i], b[(i+2)%7]*b[i]] for i in range(7)]]
    gens_labels = [[0], \
                   [i for i in range(7)], \
                   [i for i in range(7)]]
    sub_face_labels = [[(0,)], \
                       [(0,0) for i in range(7)], \
                       [((i+2)%7, (i+6)%7, i) for i in range(7)]]
    return gens,gens_labels,sub_face_labels

def LSV_quotient1():
    verts = [(p[0]._bits,p[1]._bits,p[2]._bits) for p in proj_plane(F16)]
    gens, gen_lb, subface_lb = LSV_face_generators_and_labels()
    def action(g,v):
        return tuple(coord._bits for coord in mul_mat_by_proj_vec(g, Vector(F16, vec=(F16(v[0]),F16(v[1]),F16(v[2])))))
    return Schrier_complex(verts, gens, gen_lb, subface_lb, action)


def proj_plane_flags(field):
    """
    Returns an iterable enumerating all triples (v1,v2,v3) in field^3 such that L1 = span{v1} is contained in L2 = span{v2,v3} and dim(Li) = i.
    The matrices [v1] and [v2,v3] are in reduced echelon form. 
    Each pair (L1,L2) is produced exactly once.
    """
    one = field.one()
    zero = field.zero()
    for a in field.enum():
        for b in field.enum():
            for c in field.enum():
                yield ((one,a,b), (one,zero,b-a*c), (zero,one,c))
            yield ((one,a,b), (one,a,zero), (zero,zero,one))
    for a in field.enum():
        for b in field.enum():
            yield ((zero,one,a), (one,zero,b), (zero,one,a))
        yield ((zero,one,a), (zero,one,zero), (zero,zero,one))
    for a in field.enum():
        yield ((zero,zero,one), (one,a,zero), (zero,zero,one))
    yield ((zero,zero,one), (zero,one,zero), (zero,zero,one))

def F16_flag_to_int_flag(f):
    return (f[0][0]._bits, f[0][1]._bits, f[0][2]._bits, \
            f[1][0]._bits, f[1][1]._bits, f[1][2]._bits, \
            f[2][0]._bits, f[2][1]._bits, f[2][2]._bits)

def int_flag_to_F16_flag(f):
    return(Vector(F16, vec=(F16(f[0]), F16(f[1]), F16(f[2]))), \
           Vector(F16, vec=(F16(f[3]), F16(f[4]), F16(f[5]))), \
           Vector(F16, vec=(F16(f[6]), F16(f[7]), F16(f[8]))))

def mul_mat_by_flag(M,f):
    v0 = mul_mat_by_proj_vec(M, f[0]) 
    v1 = M.apply(f[1])
    v2 = M.apply(f[2])
    # perform direct Gauss elimination on v1, v2:
    field = f[0].get_field()
    N = Matrix(field, 2, 3, mat=[[vi for vi in v] for v in [v1,v2]])
    G = N.Gauss_elim()
    return v0, G.get_row(0), G.get_row(1)

def LSV_quotient2():
    verts = [F16_flag_to_int_flag(p) for p in proj_plane_flags(F16)]
    gens, gen_lb, subface_lb = LSV_face_generators_and_labels()
    def action(g,f):
        return F16_flag_to_int_flag(mul_mat_by_flag(g, int_flag_to_F16_flag(f)))
    return Schrier_complex(verts, gens, gen_lb, subface_lb, action)

def func_from_LSVq2_to_LSVq1(face):
    return tuple(x[:3] for x in face[:-1]) + (face[-1],)

def LSV_quo_morphism_2_to_1():
    Ram1 = LSV_quotient1()
    Ram2 = LSV_quotient2()
    return DeltaSetMorphism(Ram2, func_from_LSVq2_to_LSVq1, Ram1)


def map_GLF16_Z3(g):
    det = (g[0,0]*g[1,1]*g[2,2] + g[0,1]*g[1,2]*g[2,0] + g[0,2]*g[1,0]*g[2,1]) - \
          (g[0,0]*g[1,2]*g[2,1] + g[0,1]*g[1,0]*g[2,2] + g[0,2]*g[1,1]*g[2,0])
    res = det*det
    res = res*res*det
    if res._bits == 7:
        return 1
    if res._bits == 6:
        return 2
    if res._bits == 1:
        return 0
    raise "Matrix not invertible"
    
def LSV_quotient3():
    verts = [F16_flag_to_int_flag(p)+(j,) for p in proj_plane_flags(F16) for j in range(3)]
    gens, gen_lb, subface_lb = LSV_face_generators_and_labels()
    def action(g,f):
        return F16_flag_to_int_flag(mul_mat_by_flag(g, int_flag_to_F16_flag(f[:9]))) + ((map_GLF16_Z3(g) + f[9]) % 3,)
    return Schrier_complex(verts, gens, gen_lb, subface_lb, action)

def func_from_LSVq3_to_LSVq1(face):
    return tuple(x[:3] for x in face[:-1]) + (face[-1],)

def LSV_quo_morphism_3_to_1():
    Ram1 = LSV_quotient1()
    Ram3 = LSV_quotient3()
    return DeltaSetMorphism(Ram3, func_from_LSVq3_to_LSVq1, Ram1)



