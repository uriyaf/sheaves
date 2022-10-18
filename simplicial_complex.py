
from random import randint, seed



class DeltaSet(object):
    """
    Implements a delta set, which is like a simplicial complex with ordered daces. Double edges and the like are allowed.
    The data of the delta set consists of its faces (except the empty face) and a dictionary a face f and i in {0,...,dim f} to the face obtained by removing the i-th vertex from f.
    """
    def __init__(self,faces,face_maps,empty_face):
        """
        @faces is an array of arrays. The i-th array should be the list of i-faces, consisting of immutables.
        @face_maps is an array of dictionaries. The i-th dictionary maps (f,j) with f in faces[i] and j in {0,...,i}
            to the face obtained from f by removing the i-th vertex.
        The validity of the data is not checked, but you can use validate() to check it.
        """
        self._faces = faces
        self._face_maps = face_maps
        self._face_dims = {}
        for i,layer in enumerate(faces):
            for f in layer:
                self._face_dims[f] = i
        self._face_dims[empty_face] = -1
        self._empty_face = empty_face

    def get_empty_face(self):
        return self._empty_face
        
    def count_overfaces(self,dim):
        res = {}
        for f in self.get_faces(dim):
            res[f] = 0
        for f in self.get_faces(dim+1):
            for i in range(dim+2):
                res[self._face_maps[dim+1][f,i]] += 1
        return res
            
    def get_faces(self,dim):
        """
        returns an iterator for the faces of dimension 'dim' (dim >= 0).
        """
        if 0 <= dim and dim < len(self._faces):
            for f in self._faces[dim]:
                yield f
        if dim == -1:
            yield self._empty_face # common symbol for the empty face

    def get_all_faces(self):
        yield self._empty_face
        for layer in self._faces:
            for f in layer:
                yield f

    def get_layers(self,low,high):
        """
        returns an iterator for the faces of dimension >= 'low' and smaller than 'high'.
        It is assumed that 0 <= low <= high
        """
        for layer in self._faces[low:high]:
            for f in layer:
                yield f

    def dim(self):
        return len(self._faces)-1

    def face_dim(self,f):
        return self._face_dims[f]

    def subface_codim1(self,f,i):
        return self._face_maps[self._face_dims[f]][f,i]

    def face_map(self,i):
        """
        Returns a function which is the i-th face map (i.e. remove i-th vertex from face).
        """
        return lambda f: self._face_maps[f,i]

    def subface(self,f,remove_verts):
        """
        Returns the subface of @f obtained by removing the vertices in @remove_verts.
        For example, if @remove_verts is [0,0], then the 0th vertex of @f is removed, and the the 0th vertex of the result.
        """
        res = f
        dim = self._face_dims[f]
        # remove_verts = sorted(remove_verts, reverse=True) # best not to do this
        for i in remove_verts:
            res = self._face_maps[dim][res,i]
            dim -= 1
        return res

    def validate(self):
        dim = 1
        for f in self._faces[0]:
            if self._face_maps[0][f,0] != self._empty_face:
                self._validation_failure = (f,0)
                return False
        for layer in self._faces[1:]:
            for f in layer:
                for i in range(dim+1):
                    fi = self._face_maps[dim][f,i]
                    if fi not in self._faces[dim-1]:
                        self._validation_failure = (f,i)
                        return False
                    for j in range(i):
                        fj = self._face_maps[dim][f,j]
                        fij = self._face_maps[dim-1][fi,j]
                        fji = self._face_maps[dim-1][fj,i-1]
                        if fij != fji:
                            self._validation_failure = (f,j,i)
                            return False
            dim += 1
        return True

    def __repr__(self):
        return "Simplicial complex:\n" + "\n".join(["%d-faces: %d" % (i,len(layer)) for i,layer in enumerate(self._faces)])

    def rand_face(self,dim):
        """
        Returns a random face, chosen uniformly among the faces for the given dimension.
        It is assumed that @dim >= 0.
        """
        return self._faces[dim][randint(0,len(self._faces[dim])-1)]

##    def _add_free_triangle(self,edge,v0,e01,e02,t):
##        # do not use
##        e12 = edge
##        v2 = self._face_maps[1][e12,0]
##        v1 = self._face_maps[1][e12,1]
##        empty = self._face_maps[0][v1,0]
##        self._faces[0].append(v0)
##        self._faces[1].append(e01)
##        self._faces[1].append(e02)
##        self._faces[2].append(t)
##        self._face_dims[v0] = 0
##        self._face_dims[e01] = 1
##        self._face_dims[e02] = 1
##        self._face_dims[t] = 2
##        self._face_maps[0][v0,0] = empty
##        self._face_maps[1][e01,0] = v1
##        self._face_maps[1][e01,1] = v0
##        self._face_maps[1][e02,0] = v2
##        self._face_maps[1][e02,1] = v0
##        self._face_maps[2][t,0] = e12
##        self._face_maps[2][t,1] = e02
##        self._face_maps[2][t,2] = e01


class DeltaSetMorphism(object):
    """
    Implements a morphism between DeltaSets
    """
    def __init__(self,source,func,target):
        """
        @source and @target and DeltaSets
        @func is a function from faces of @source to faces of @target.
        Use the method validate() to vadiate the morphism.
        """
        self.source = source
        self.func = func
        self.target = target

    def __call__(self,face):
        """
        Apply morphism to a face.
        """
        return self.func(face)

    def validate(self):
        for i in range(self.source.dim()+1):
            for j in range(i+1):
                for f in self.source.get_faces(i):
                    if self.func(self.source._face_maps[i][f,j]) != self.target._face_maps[i][self.func(f),j]:
                        return False
        return True
    
def rename_faces_as_integers(sc):
    """
    Takes a simplicial complex (i.e. a DeltaSet) and returns a pair consisting of:
    (1) an isomorphic simplicial complex (DeltaSet) whose faces are labelled by integers,
    (2) and isomorphism (DelatSetMorphism) from @sc to the new simplicial complex
    """
    rename = {sc.get_empty_face():-1}
    name = 0
    faces = []
    for i in range(sc.dim()+1):
        layer = []
        for f in sc.get_faces(i):
            rename[f] = name
            layer.append(name)
            name += 1
        faces.append(layer)
    face_maps = [{}]
    for f in sc.get_faces(0):
        face_maps[0][rename[f],0]=-1
    for i in range(1,sc.dim()+1):
        layer = {}
        for f in sc.get_faces(i):
            for j in range(i+1):
                layer[rename[f],j] = rename[sc._face_maps[i][f,j]]
        face_maps.append(layer)
    target = DeltaSet(faces, face_maps, -1)
    mor = DeltaSetMorphism(sc, lambda x: rename[x], target)
    return target, mor


# the following functions generate some examples of simplicial complexes.

def Schrier_complex(vertices, face_generators, face_labels, sub_face_labels, action, sort_verts=True):
    """
    A general function for constructing the Schrier complex from a list of "genertors" acting on a
    list of "labels". Takes double faces into account (which is why @sub_face_labels are needed).
    The input is as follows:
    @vertices is a list of hashables that will be the vertices of the simplicial complex.
    @action is a function taking a vertex and a generator and outputing another vertex.
    @face_generators is a list of lists. Each entry of the i-th list must be list [g0,g1,...,gi] of generators
        indicating that the i-face (g0*v, g1*v, ..., gi*v) must be included in the complex.
    @face_labels is a list of lists. The i-th list must consist of distinct hashables. The j-th hashable is the
        label of the j-th face indictor in @face_generators[i]. For example, if @face_generators[i][j] = [g0,g1,...,gi],
        then the face corresponding to this list of generators and a vertex v will be named
        (g0, ... , gi, @face_labels[i][j]); the label is important if the Schrier complex has
        double faces. 
    @sub_face_labels is a list of lists of lists. The list at the (i,j)-place must consist of i+1 face labels
        l=[l[0],...,l[i]] taken from amount @face_labels[i-1], with l[k] indicating the label of the face obtained
        by removing the k-th vertex (starting from 0) at the @face_generators[i][j] face. (For example, if
        @face_generators[i][j]=[g0,g1,...,gi], then the label of the face
        (g0*v, g1*v, ..., g(k-1)*v, ..., g(k+1)*v, ...., gi*v) will be l[k].)
    @sort_verts if set to true, then when naming the faces of the simplicial complex, the list (g0*v, g1*v, ..., gi*v)
        will be sorted. (This may jumble the labels.)
    """
    faces = []
    face_maps = []
    empty_face = tuple(sub_face_labels[0][0])
    for i in range(len(face_generators)):
        ifaces = []
        iface_map = {}
        for v in vertices:
            for fg, label, subface_lb in zip(face_generators[i], face_labels[i], sub_face_labels[i]):
                orig_verts = [action(g,v) for g in fg]
                if sort_verts:
                    face = tuple(sorted(orig_verts)+[label])
                    ifaces.append(face)
                    for j in range(i+1):
                        subface = tuple(sorted(orig_verts[:j]+orig_verts[j+1:])+[subface_lb[j]])
                        iface_map[face,face.index(orig_verts[j])] = subface
                else:
                    face = tuple(orig_verts+[label])
                    ifaces.append(face)
                    for j in range(i+1):
                        subface = tuple(orig_verts[:j]+orig_verts[j+1:]+[subface_lb[j]])
                        iface_map[face,j] = subface
        faces.append(ifaces)
        face_maps.append(iface_map)
    return DeltaSet(faces, face_maps, empty_face)


def torus_2dim_gen_and_rels():
    gens = [[((0,0),)], \
            [((0,0), x) for x in [(0,1),(1,0),(1,1)]], \
            [[(0,0),(0,1),(1,1)], [(0,0),(1,0),(1,1)]]]
    gens_labels = [[0], \
                   [0,1,2], \
                   [0,1]]
    sub_face_labels = [[(0,)], \
                       [(0,0), (0,0), (0,0)], \
                       [(1,2,0), (0,2,1)]]
    return gens,gens_labels,sub_face_labels

def torus_2dim(grid=(3,3)):
    verts = [(i,j) for i in range(grid[0]) for j in range(grid[1])]
    gens, gen_lb, subface_lb = torus_2dim_gen_and_rels()
    def action(g,f):
        return ((g[0]+f[0]) % grid[0], (g[1]+f[1]) % grid[1])
    return Schrier_complex(verts, gens, gen_lb, subface_lb, action, False)

def torus_2dim_morphism(grid=(3,3), covering_factors = (3,3)):
    T = torus_2dim(grid)
    S = torus_2dim((grid[0]*covering_factors[0], grid[1]*covering_factors[1]))
    def S_to_T(f):
        return tuple((x % grid[0], y % grid[0]) for x,y in f[:-1]) + (f[-1],)
    return DeltaSetMorphism(S, S_to_T, T)


def sum_tuple(L1,L2):
    return tuple(x+y for x,y in zip(L1,L2))

def sub_tuple(L1,L2):
    return tuple(x-y for x,y in zip(L1,L2))

def all_01_tuples(length):
    """
    Returns all tuples of length 'length', ordered lexicographically.
    """
    if length <= 0:
        return [()]
    T = all_01_tuples(length-1)
    return [(0,)+t for t in T] + [(1,)+t for t in T]


def torus_ndim_gens_and_rels(dim, skeleton=None):
    """
    Constructs a triangulation of the n-torus with n=len(dims).
    If dims = (n1,n2,...) then the triangulation has vertices (x1,x2,...) with xi ranging on {0,...,ni-1}.
    Only the skeleton up to dimension 'skeleton' (including) is constructed.
    """
    n = dim
    if skeleton is None:
        skeleton = n
    if dim<0 or skeleton<0:
        raise "Dimension should be at least 0."
    leaps = all_01_tuples(n)[1:]
    vertex_gens = [[tuple(0 for d in range(n))]]
    all_gens = [vertex_gens]
    for i in range(1,skeleton+1):
        new_gens = []
        for l in all_gens[-1]:
            for t in leaps:
                s = sum_tuple(l[-1],t)
                if 2 not in s:
                    new_gens.append(l+[s])
        all_gens.append(new_gens)
    labels = [[j for j in range(len(gens))] for gens in all_gens]
    subface_labels = [[(0,)]]
    for i in range(1,skeleton+1):
        new_labels = []
        for leaps in all_gens[i]:
            sfl = []
            for j in range(len(leaps)):
                l = leaps[:j] + leaps[j+1:]
                lead = l[0]
                l = [sub_tuple(li,lead) for li in l]
                ind = all_gens[i-1].index(l)
                sfl.append(ind)
            new_labels.append(tuple(sfl))
        subface_labels.append(new_labels)
    return all_gens, labels, subface_labels

def generate_grid(dims):
    if len(dims)==0:
        return [()]
    T = generate_grid(dims[1:])
    return sum(([(i,)+t for t in T] for i in range(dims[0])), [])
            
def torus_ndim(grid,skeleton=None):
    """
    Returns a simplicial complex which is a triangulation of an n-dimensional-torus.
    Here, n = len(@grid).
    @grid is an array of positive integers, specifiying the dimensions of the grid,
        e.g., @grid=[3,2,2] will give a torus triangulation which refines a 3-by-2-by-2
        box grid.
    """
    verts = generate_grid(grid)
    gens, gen_lb, subface_lb = torus_ndim_gens_and_rels(len(grid),skeleton)
    def action(g,f):
        return tuple((g[i]+f[i]) % grid[i] for i in range(len(grid)))
    return Schrier_complex(verts, gens, gen_lb, subface_lb, action, False)

def torus_ndim_morphism(grid=(3,3,3),factors=(3,3,3),skeleton=None):
    T = torus_ndim(grid,skeleton)
    S = torus_ndim(tuple(i*j for i,j in zip(grid,factors)), skeleton)
    def S_to_T(f):
        return tuple(tuple(x % g for x,g in zip(fi,grid)) for fi in f[:-1]) + (f[-1],)
    return DeltaSetMorphism(S, S_to_T, T)

