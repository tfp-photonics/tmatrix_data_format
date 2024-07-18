"""A small tool to parse different types of mesh files and convert them.

This tool provides utilites to parse mesh files. The mesh file can then be written to
a different format. The goal is to reproduce the content of the original mesh as
faithful as possible. However, an exact one to one correspondance is not always
possible. Currently, COMSOL files (.mphtxt), JCMsuite files (.jcm), and gmsh files
(.msh) are supported.
"""

import operator

import argparse as _argparse
import io as _io
import os as _os
import random as _random
import struct as _struct
from abc import ABCMeta as _ABCMeta
from collections import abc as _abc
from importlib import metadata as _metadata
import numpy as np

# __version__ = _metadata.version("meshparser")
_MIRROR = 2**28
_TRANSPARENT = 2**29
_PERIODIC = 2**30
_MIRROR_JCM = -29_856_126


class Node(_abc.Sequence):
    """Node of a mesh.

    A node is a single point in 2D or 3D space. It is associated with a tag. The node is
    defined either as a point in 2D space with the `x` and `y` coordinates or in 3D
    space with `x`, `y`, and `z` coordinates. A node is the most basic constituent of a
    mesh. Nodes behave similar to tuples of length two or three: they are iterable,
    immutable, and hashable. A comparison also compares tags.
    """

    def __init__(self, coordinate, *, tag):
        """Initialize node.

        Args:
            coordinate (Iterable): coordinates defining the x, y(, and z) component
            tag (int): identifier
        """
        self._tag = int(tag)
        self._x, self._y, *z = map(float, coordinate)
        self._sdim = 2 + len(z)
        if self._sdim == 2:
            self._z = 0.0
        elif self._sdim == 3:
            self._z = z[0]
        else:
            raise ValueError(f"space dimension is '{self._sdim}' but expected 2 or 3")

    @property
    def x(self):
        """float: X coordinate."""
        return self._x

    @property
    def y(self):
        """float: Y coordinate."""
        return self._y

    @property
    def z(self):
        """float: Z coordinate."""
        return self._z

    @property
    def tag(self):
        """int: Identifier."""
        return self._tag

    @property
    def sdim(self):
        """int: Spatial dimension."""
        return self._sdim

    @property
    def in_3d(self):
        """tuple: Coordinates embedded in 3D space."""
        return self._x, self._y, self._z

    def __repr__(self):
        return f"Node({self[:]}, tag={self.tag})"

    def __len__(self):
        return self._sdim

    def __getitem__(self, key):
        """ """
        if self.sdim == 2:
            return (self.x, self.y)[key]
        return (self.x, self.y, self.z)[key]

    def __hash__(self):
        return hash(self[:] + (self.tag,))

    def __eq__(self, other):
        return self.tag == other.tag and self[:] == other[:]


class _ElementMeta(_ABCMeta):
    """Metaclass for :class:`Element`.

    Require the attributes `dim` and `n_nodes`. If not specified, derive the COMSOL and
    JCMsuite names from the class name. The metaclass has to inherit from `ABCMeta` to
    keep the abstract base class `Sequence` working.
    """

    def __new__(cls, clsname, bases, attrs, /, **kwargs):
        for attr in ("dim", "n_nodes"):
            if attr not in attrs:
                TypeError(f"can't create class {clsname} without attribute {attr}.")
        dim = attrs["dim"]
        attrs.setdefault("comsol_name", "3 " + clsname[:3].lower())
        jcm_name = attrs.setdefault("jcm_name", clsname)
        attrs.setdefault("jcm_short_name", "#" + jcm_name[:2])
        attrs["jcm_info"] = {
            dim: f"<I>N{jcm_name}s",
            dim + 1: f"<I>NBoundary{jcm_name}s",
        }
        pre_bulk = f"# {jcm_name}s ###################\n# {jcm_name}s (p"
        pre_boundary = f"# {jcm_name} (p"
        middle = ", p".join(map(str, range(1, attrs["n_nodes"] + 1)))
        post_bulk = ", id)\n"
        post_boundary = ", boundary id, topology indicator, (identic points))\n"
        attrs["jcm_header"] = {
            dim: pre_bulk + middle + post_bulk,
            dim + 1: pre_boundary + middle + post_boundary,
        }
        return super().__new__(cls, clsname, bases, attrs, **kwargs)


class Element(_abc.Sequence, metaclass=_ElementMeta):
    """Mesh element.

    Elements are a collection of :class:`Node` that define a specific shape. The element
    itselfs defines the shape and it is assumed that the nodes are set accordingly. The
    dimension of the element is the dimension of the shape not the spatial dimension.
    Elements are iterable, immutable, and hashable. A comparison also compares tags.
    each class should specify attributes that describe how they are represented in
    different mesh files.
    """

    dim = -1
    n_nodes = 0

    def __init__(self, nodes, *, tag):
        """Initialize element.

        Args:
            nodes (Iterable): nodes of the element
            tag (int): identifier
        """
        self._sdim = nodes[0].sdim
        if self.sdim < self.dim:
            raise ValueError("node dimension must be larger than element dimension")
        for node in nodes:
            if node.sdim != self.sdim:
                raise ValueError("element contains nodes of mixed dimension")
        self._nodes = tuple(nodes)
        self._tag = int(tag)
        self.jcm_name = type(self).__name__

    @property
    def tag(self):
        """int: Identifier."""
        return self._tag

    @property
    def sdim(self):
        """int: Spatial dimension."""
        return self._sdim

    def __repr__(self):
        return f"{type(self).__name__}({self._nodes}, tag={self.tag})"

    def __len__(self):
        return self.n_nodes

    def __getitem__(self, key):
        return self._nodes[key]

    def __hash__(self):
        return hash((self.__class__, self.tag) + self._nodes)

    def __eq__(self, other):
        return (
            self.__class__ == other.__class__
            and self.tag == other.tag
            and all(a == b for a, b in zip(self, other))
        )


class Vertex(Element):
    """Actually just a node.

    The vertex gives another name to a node. As such it is a 0D mesh element. It exists
    in COMSOL meshes. In gmsh meshes it exists somewhat implicity as `Point`.
    """

    dim = 0
    n_nodes = 1
    # It gets an gmsh type of None because it exists in gmsh only as `Point`, which is
    # an Entity (i.e. a Domain) but not as Element
    gmsh_type = None
    comsol_name = "3 vtx"


class Edge(Element):
    """Two nodes connected.

    An edge is a 1D mesh element. It exists in JCMsuite only as a boundary element of
    2D meshes.
    """

    dim = 1
    n_nodes = 2
    gmsh_type = 1
    jcm_short_name = "#E"


class Triangle(Element):
    """Three nodes that are not on a straight line.

    Triangles are the central element for 2D meshes.
    """

    dim = 2
    n_nodes = 3
    gmsh_type = 2


class Quadrilateral(Element):
    """Four nodes in a plane.

    The order of the nodes is different in different programs. In JCMsuite and gmsh they
    are in clockwise order, in COMSOL the last two nodes are exchanged.
    """

    dim = 2
    n_nodes = 4
    gmsh_type = 3
    comsol_name = "4 quad"
    jcm_short_name = "#Q"


class Tetrahedron(Element):
    """Four nodes that are not in a plane."""

    dim = 3
    n_nodes = 4
    gmsh_type = 4


class Pyramid(Element):
    """Five nodes forming a quadrilateral and four triangles.

    Like for the :class:`Quadrilateral`, COMSOL uses a different node order for the
    base. The last node is always the pyramid top.
    """

    dim = 3
    n_nodes = 5
    gmsh_type = 7


class Prism(Element):
    """Six nodes forming two triangles separated by some distance."""

    dim = 3
    n_nodes = 6
    gmsh_type = 6
    comsol_name = "5 prism"


class Hexahedron(Element):
    """Eight nodes forming two quadrilaterals separated by some distance.

    Like for the :class:`Quadrilateral`, COMSOL defines them with a different node
    order. JCMsuite, inconsistently, does so as well.
    """

    dim = 3
    n_nodes = 8
    gmsh_type = 5
    jcm_name = "Brick"
    jcm_short_name = "#B"


class Domain(_abc.Set):
    """Domain of a mesh.

    A domain is the highest constituent of a mesh. It collects multiple elements of the
    same dimension that form a common shape or physical object together. It is
    essentially an immutable (and, thus, hashable) set. The order of the elements
    doesn't matter for the set operations, however, the elements are sorted by their tag
    upon intialization by their tag, to provide a fixed order when writing the mesh to a
    file. A domain can be empty, then its dimension and spatial dimension have to be
    given. In gmsh domains are named `Entity`.
    """

    def __init__(self, elements, *, tag, dim=None, sdim=None):
        """Initialize domain.

        Args:
            elements (Iterable): constituents of the domain
            tag (int): identifier
            dim (int, optional): dimension
            sdim (int, optional): spation dimension
        """
        if len(elements) != len(set(elements)):
            raise ValueError("duplicate elements")
        self._elements = sorted(elements, key=lambda item: item.tag)
        dims = {e.dim for e in self}
        if dim is not None:
            dims.add(dim)
        sdims = {e.sdim for e in self}
        if sdim is not None:
            sdims.add(sdim)
        if len(dims) != 1 or len(sdims) != 1:
            raise ValueError("invalid dimensions")
        self._dim = next(iter(dims))
        self._sdim = next(iter(sdims))
        self._tag = int(tag)

    @property
    def tag(self):
        """int: Identifier."""
        return self._tag

    @property
    def dim(self):
        """int: Dimension of the domain."""
        return self._dim

    @property
    def sdim(self):
        """int: Spatial dimension."""
        return self._sdim

    def __repr__(self):
        return "Domain({...}, " f"tag={self.tag}, dim={self.dim}, sdim={self.sdim})"

    def __iter__(self):
        return iter(self._elements)

    def __len__(self):
        return len(self._elements)

    def __contains__(self, item):
        return item in self._elements

    def __hash__(self):
        return hash((self._hash(), self.tag))

    def __le__(self, other):
        return self.tag == other.tag and super().__le__(other)

    def __ge__(self, other):
        return self.tag == other.tag and super().__le__(other)


class MeshError(Exception):
    """Error for :class:`Mesh`."""


class Mesh:
    """Mesh definition.

    A mesh is the collection of one or multiple domains, their elements and nodes. The
    mesh must be 'closed', so all nodes must belong to an element and the nodes of all
    elements must be included. Similarly, this applies to elements and domains. Changes
    of the attributes are not forbidden so changes are at the users discretion. However,
    the consistency of the attributes is as well.

    The domains, elements, and nodes are stored as keys in dictionaries with an integer
    tag (or None) as value. This tag can be subject to change. However, the 'original'
    tag can be restored from the constituents itself.

    Attributes:
        domains (dict): :class:`Domain`, tag (int) pairs
        elements (dict): :class:`Element`, tag (int) pairs
        nodes (dict): :class:`Node`, tag (int) pairs
        dims (set): dimensions of the domains and elements
        sdim (int): spatial dimension
        origin (str): type of mesh from which it was loaded ('comsol', 'jcm', 'gmsh')
        info (dict): additional information
    """

    def __init__(
        self, domains=None, elements=None, nodes=None, *, origin=None, info=None
    ):
        """Intialize mesh.

        Args:
            domains (dict): :class:`Domain`, tag (int) pairs
            elements (dict): :class:`Element`, tag (int) pairs
            nodes (dict): :class:`Node`, tag (int) pairs
            origin (str): type of mesh from which it was loaded ('comsol', 'jcm',
                'gmsh')
            info (dict): additional information
        """
        if domains is None:
            domains = {}
        if elements is None:
            elements = {e: e.tag for domain in domains for e in domain}
        if nodes is None:
            nodes = {n: n.tag for element in elements for n in element}
        self.sdim = next(iter(nodes)).sdim
        self.dims = {domain.dim for domain in domains}
        self.nodes = dict(nodes)
        self.elements = dict(elements)
        self.domains = dict(domains)
        self._sdim = self.sdim  # if sdim is 3 and all z values are zero, this is 2
        self.check()
        self.origin = origin
        self.info = {} if info is None else info

    def __repr__(self):
        return (
            "Mesh({...}, {...}, {...}, "
            f"dims={self.dims}, sdim={self.sdim}, origin='{self.origin}')"
        )

    def check(self):
        """Check the mesh for consistency."""
        elements = {element for domain in self.domains for element in domain}
        if elements != self.elements.keys():
            raise MeshError("'domains' and 'elements' incompatible")
        nodes = {node for element in self.elements for node in element}
        if nodes != self.nodes.keys():
            raise MeshError("'elements' and 'nodes' incompatible")
        for node in self.nodes:
            if node.sdim != self.sdim:
                raise MeshError("mesh contains nodes of mixed dimension")
        for domain in self.domains:
            if domain.dim not in self.dims:
                raise MeshError(f"'dims' is missing dimension '{domain.dim}'")
        if self.sdim == 3 and max(self.dims) < 3 and all(n[2] == 0 for n in self.nodes):
            self._sdim = 2

    def elements_by_type(self):
        """Return the elements sorted by their type.

        Returns:
            dict
        """
        ret = {}
        for element, tag in self.elements.items():
            ret.setdefault(element.__class__, {})[element] = tag
        return {e: ret[e] for e in Element.__subclasses__() if e in ret}

    def tag_restore(self):
        """Update dictionaries with the tags of each constituent."""
        for domain in self.domains:
            self.domains[domain] = domain.tag
        for element in self.elements:
            self.elements[element] = element.tag
        for node in self.nodes:
            self.nodes[node] = node.tag

    def tag_comsol(self):
        """Tag domains to be compatible with COMSOL."""
        if self.origin == "comsol":
            self.tag_restore()
            return
        self.nodes = {n: i for i, n in enumerate(self.nodes)}
        domain_tags = [0] * (self._sdim + 1)
        domain_tags[self._sdim] = 1
        jcm_boundary = self.info.get("jcm_boundary", {})
        for domain in self.domains:
            tag = domain_tags[domain.dim]
            domain_tags[domain.dim] += 1
            mode = set(self.info.get("phystags", {}).get(domain, [])) & {
                _MIRROR,
                _PERIODIC,
                _TRANSPARENT,
            }
            if len(mode) > 1:
                raise MeshError("unrecognized special mode in 'phystags'")
            if len(mode) == 1:
                mode = next(iter(mode))
            else:
                mode = None
            for element in domain:
                top = int(jcm_boundary.get(element, [0])[0])
                if mode == _MIRROR or element.tag == _MIRROR_JCM:
                    self.elements[element] = tag + _MIRROR
                elif mode == _PERIODIC or top == -1:
                    self.elements[element] = tag + _PERIODIC
                elif mode == _TRANSPARENT or top == -2:
                    self.elements[element] = tag + _TRANSPARENT
                else:
                    self.elements[element] = tag

    def tag_jcm(self):
        """Tag domains to be compatible with JCMsuite."""
        if self.origin == "jcm":
            self.tag_restore()
            return
        self.nodes = {n: i + 1 for i, n in enumerate(self.nodes)}
        domain_tags = [1] * (self._sdim + 1)
        self.info["jcm_boundary"] = {}  # delete previous versions
        for domain in self.domains:
            tag = domain_tags[domain.dim]
            domain_tags[domain.dim] += 1
            mode = set(self.info.get("phystags", {}).get(domain, [])) & {
                _MIRROR,
                _PERIODIC,
                _TRANSPARENT,
            }
            if len(mode) > 1:
                raise MeshError("unrecognized special mode in 'phystags'")
            if len(mode) == 1:
                mode = next(iter(mode))
            else:
                mode = None
            for element in domain:
                if mode == _MIRROR or element.tag & _MIRROR:
                    self.elements[element] = _MIRROR_JCM
                elif mode == _PERIODIC or element.tag & _PERIODIC:
                    # Actually, here we would need to define the jcm_boundary -1 and the
                    # nodes on the opposite side. However, we ignore this info for now.
                    self.elements[element] = tag
                elif mode == _TRANSPARENT or element.tag & _TRANSPARENT:
                    self.elements[element] = tag
                    self.info["jcm_boundary"][element] = ["-2"]
                else:
                    self.elements[element] = tag

    def tag_gmsh(self):
        """Tag domains to be compatible with gmsh."""
        if self.origin == "gmsh":
            self.tag_restore()
            return
        self.nodes = {n: i + 1 for i, n in enumerate(self.nodes)}
        self.elements = {n: i + 1 for i, n in enumerate(self.elements)}
        self.domains = {n: i + 1 for i, n in enumerate(self.domains)}
        jcm_boundary = self.info.get("jcm_boundary", {})
        self.info["phystags"] = {}  # delete previous versions
        for domain in self.domains:
            if len(domain) == 0:
                continue
            element = next(iter(domain))
            top = jcm_boundary.get(element, [0])[0]
            if element.tag == _MIRROR_JCM or element.tag & _MIRROR:
                self.info["phystags"][domain] = [_MIRROR]
            elif top == -1 or element.tag & _PERIODIC:
                self.info["phystags"][domain] = [_PERIODIC]
            elif top == -2 or element.tag & _TRANSPARENT:
                self.info["phystags"][domain] = [_TRANSPARENT]

    @classmethod
    def read(cls, name):
        """Read file.

        The file with the given name is read. The type of the file is determined by its
        extentsion.

        Args:
            name (str): file name
        """
        funcs = {
            ".jcm": cls.read_jcm,
            ".mphtxt": cls.read_comsol,
            ".msh": cls.read_gmsh,
            ".msh4": cls.read_gmsh,
        }
        func = funcs[_os.path.splitext(name)[-1]]
        with open(name, mode="rb") as fobj:
            return func(fobj)

    def write(self, name, mode="x", *, skip_1d=False):
        """Write the mesh to a file.

        The mesh is written to a file with the given name. The mode can be chosen
        according to the documentation of `open()`. The legacy JCMsuite mesh files can
        also be written (applies to 2D meshes only).

        Args:
            name (str): file name
            mode (str, default 'x'): writing mode (overwrite 'w', append 'a', write but
                fail when already exiting 'x')
            skip_1d (bool, default False): Skip writing 1D elements for 2D meshes in
                COMSOL, because they can cause loading errors there. They can be kept,
                if all 2D domains are bounded by 1D elements.
        """
        funcs = {
            ".jcm": self.write_jcm,
            ".mphtxt": self.write_comsol,
            ".msh": self.write_gmsh,
            ".msh4": self.write_gmsh,
        }
        func = funcs[_os.path.splitext(name)[-1]]
        kwargs = {}
        if "b" in mode and func != self.write_jcm:
            raise ValueError("invalid mode")
        if skip_1d and func == self.write_comsol:
            kwargs["skip_1d"] = True
        with open(name, mode) as fobj:
            func(fobj, **kwargs)

    @classmethod
    def read_jcm(cls, fobj):
        """Read a JCMsuite grid file (.jcm).

        Given a file object of the jcm file, it is read and the corresponding mesh is
        created.

        Args:
            fobj (file object): file to read
        """
        header = _read_header_jcm(fobj)
        dim = int(header["<I>SpaceDim"])
        if dim == int(header["<I>ManifoldDim"]) == 2:
            element_types = (Triangle, Quadrilateral, Edge)
        elif dim == int(header["<I>ManifoldDim"]) == 3:
            element_types = (
                Tetrahedron,
                Pyramid,
                Prism,
                Hexahedron,
                Triangle,
                Quadrilateral,
            )
        else:
            raise ValueError("invalid dimensions")

        domains_bulk = {}
        domains_boundary = {}
        nodes = {}
        elements = {}
        boundary_info = {}

        if header["__MODE__"] == "BINARY0":
            if isinstance(fobj, _io.TextIOBase):
                raise ValueError("binary file must be opened in binary mode")
            if fobj.read(1) != b"0":
                raise ValueError("invalid binary file")
            binary = True
        elif header["__MODE__"] == "TEXT":
            fobj = _iter_wo_comments(fobj)
            binary = False
        else:
            raise ValueError("invalid '__MODE__'")
        for _ in range(int(header["<I>NPoints"])):
            tag, *values = _read_jcm(fobj, f"i{dim}d", binary)
            nodes[tag] = Node(values, tag=tag)
        for ecls in element_types:
            for _ in range(int(header[ecls.jcm_info[dim]])):
                *values, tag = _read_jcm(fobj, f"{ecls.n_nodes + 1}i", binary)
                element = ecls(
                    _jcm_node_order([nodes[i] for i in values], ecls),
                    tag=tag,
                )
                elements[element] = tag
                if ecls.dim == dim:
                    domains_bulk.setdefault(tag, []).append(element)
                else:
                    top = _read_jcm(fobj, "i", binary)
                    if top[0] == -1:
                        top = top + _read_jcm(fobj, f"{ecls.n_nodes}i", binary)
                    boundary_info[element] = top
                    domains_boundary.setdefault((tag, top[0]), []).append(element)
        domains = {Domain(val, tag=tag): tag for tag, val in domains_bulk.items()}
        domains.update(
            {Domain(val, tag=tag): tag for (tag, _), val in domains_boundary.items()}
        )
        return cls(
            domains,
            elements,
            {v: k for k, v in nodes.items()},
            origin="jcm",
            info={"jcm_boundary": boundary_info, "jcm_header": header},
        )

    @classmethod
    def read_comsol(cls, fobj):
        """Read a COMSOL mesh file (.mphtxt).

        Given a file object of the mphtxt file, it is read and the corresponding mesh is
        created.

        Args:
            fobj (file object): file to read
        """
        # Domains for 0d, 1d, 2d, and 3d
        domains = {}

        element_types = {e.comsol_name: e for e in Element.__subclasses__()}

        itr = _iter_wo_comments(fobj)
        # 1. version
        if next(itr) != "0 1":
            raise ValueError("invalid file version number")
        # 2. skip tags
        meshname = " \n".join(next(itr) for _ in range(int(next(itr))))
        # 3. skip types
        for _ in range(int(next(itr))):
            next(itr)
        # 4. object number?
        if next(itr) != "0 0 1":
            raise ValueError("invalid file object number")
        # 5. class
        if next(itr) != "4 Mesh":
            raise ValueError("invalid class")
        # 6. another version
        if next(itr) != "4":
            raise ValueError("invalid version")
        # 7. dimension
        sdim = int(next(itr))
        # 8. nodes
        consume = int(next(itr))
        nodes = {}
        tag_min = int(next(itr))
        for i in range(consume):
            nodes[tag_min + i] = Node(map(float, next(itr).split()), tag=tag_min + i)
        # 9. Elements
        elements = {}
        for _ in range(int(next(itr))):
            mode = next(itr)
            ecls = element_types[mode]
            dim = ecls.dim
            if ecls.n_nodes != int(next(itr)):
                raise ValueError(f"error reading '{mode}'")
            new_elements = []
            for __ in range(int(next(itr))):
                values = next(itr).split()
                new_elements.append(
                    _comsol_node_order([nodes[int(i)] for i in values], ecls)
                )
            for i in range(int(next(itr))):
                tag = int(next(itr))
                element = ecls(new_elements[i], tag=tag)
                elements[element] = tag
                domains.setdefault((dim, tag), []).append(element)
        domains = {
            Domain(val, tag=tag, dim=dim, sdim=sdim): tag
            for (dim, tag), val in domains.items()
        }
        return cls(
            domains,
            elements,
            {v: k for k, v in nodes.items()},
            origin="comsol",
            info={"meshname": meshname},
        )

    @classmethod
    def read_gmsh(cls, fobj):
        """Read a gmsh mesh file (.msh).

        Given a file object of the gmsh file, it is read and the corresponding mesh is
        created.

        Args:
            fobj (file object): file to read
        """
        itr = _gmsh_iter(fobj)
        element_types = {e.gmsh_type: e for e in Element.__subclasses__()}
        domains = {}
        elements = {}
        nodes = {}
        nodes0d = {}
        physnames = {}
        phystags = {}
        boundingpoints = {}
        for line in itr:
            if line == "$MeshFormat":
                if next(itr) != "4.1 0 8":
                    raise ValueError("unsupported version")
                if next(itr) != "$EndMeshFormat":
                    raise ValueError("invalid section 'MeshFormat'")
            elif line == "$PhysicalNames":
                for i in range(int(next(itr))):
                    dim, tag, *name = next(itr).split(" ")
                    physnames[(int(dim), int(tag))] = " ".join(name)
                if next(itr) != "$EndPhysicalNames":
                    raise ValueError("invalid section 'PhysicalNames'")
            elif line == "$Entities":
                blocks = map(int, next(itr).split())
                for dim, nblocks in enumerate(blocks):
                    for i in range(nblocks):
                        values = next(itr).split()
                        tag = int(values[0])
                        start = 7
                        if dim == 0:
                            start = 4
                            # gmsh add their "Points" only here without a node, so we
                            # store them with negative tags and append them in the end
                            # to the other nodes
                            nodes0d[-tag] = Node(map(float, values[1:4]), tag=-tag)
                            element = Vertex([nodes0d[-tag]], tag=tag)
                            elements[element] = tag
                            domains[(0, tag)] = [element]
                        inc = int(values[start])
                        start += 1
                        phystags[(dim, tag)] = list(
                            map(int, values[start : start + inc])
                        )
                        if dim == 0:
                            continue
                        start += inc
                        inc = int(values[start])
                        start += 1
                        boundingpoints[(dim, tag)] = list(
                            map(int, values[start : start + inc])
                        )
                if next(itr) != "$EndEntities":
                    raise ValueError("invalid section 'Entitites'")
            elif line == "$Nodes":
                nblocks = int(next(itr).split()[0])
                for _ in range(nblocks):
                    dim, domaintag, _, ntags = map(int, next(itr).split())
                    tags = [int(next(itr)) for _ in range(ntags)]
                    for tag in tags:
                        nodes[tag] = Node(map(float, next(itr).split()), tag=tag)
                    domains.setdefault((dim, domaintag), [])
                if next(itr) != "$EndNodes":
                    raise ValueError("invalid section 'Nodes'")
            elif line == "$Elements":
                nblocks = int(next(itr).split()[0])
                for _ in range(nblocks):
                    dim, domaintag, element_type, ntags = map(int, next(itr).split())
                    ecls = element_types[element_type]
                    domain = domains.setdefault((dim, domaintag), [])
                    for _ in range(ntags):
                        tag, *nodetags = map(int, next(itr).split())
                        element = ecls([nodes[i] for i in nodetags], tag=tag)
                        elements[element] = tag
                        domain.append(element)
                if next(itr) != "$EndElements":
                    raise ValueError("invalid section 'Elements'")
        nodes.update(nodes0d)
        sdim = next(iter(nodes.values())).sdim
        domains = {
            (dim, tag): Domain(val, tag=tag, dim=dim, sdim=sdim)
            for (dim, tag), val in domains.items()
        }
        phystags = {domains[k]: v for k, v in phystags.items()}
        boundingpoints = {domains[k]: v for k, v in boundingpoints.items()}
        domains = {v: k for k, v in domains.items()}
        info = {
            "phystags": phystags,
            "boundingpoints": boundingpoints,
        }
        if physnames:
            info["physnames"] = physnames
        return cls(
            domains,
            elements,
            {v: k for k, v in nodes.items()},
            origin="gmsh",
            info=info,
        )

    def write_comsol(self, fobj, *, skip_1d=False):
        """Write a COMSOL file (.mphtxt).

        Args:
            fobj (file object): file to write to
            skip_1d (bool, default False): Skip writing 1D elements for 2D meshes in
                COMSOL, because they can cause loading errors there. They can be kept,
                if all 2D domains are bounded by 1D elements.
        """
        self.tag_comsol()
        fobj.write(
            """# Created by meshparser.

# Major & minor version
0 1 
"""  # noqa: W291
        )
        meshname = self.info.get("meshname", "5 mesh1")
        fobj.write(str(len(meshname.split("\n"))))
        fobj.write(
            f""" # number of tags
# Tags
{meshname} 
1 # number of types
# Types
3 obj 

# --------- Object 0 ----------

0 0 1 
4 Mesh # class
4 # version
{self._sdim} # sdim
{len(self.nodes)} # number of mesh vertices
0 # lowest mesh vertex index

# Mesh vertex coordinates
"""  # noqa: W291
        )
        for node in self.nodes:
            fobj.write(" ".join(f"{i:.17g}" for i in node[: self._sdim]) + " \n")

        elements_by_type = self.elements_by_type()
        if self._sdim == 2 and skip_1d:
            elements_by_type.pop(Edge, None)
        fobj.write(
            f"""
{len(elements_by_type)} # number of element types
"""
        )
        for i, (ecls, elements) in enumerate(elements_by_type.items()):
            fobj.write(
                f"""
# Type #{i}

{ecls.comsol_name} # type name


{ecls.n_nodes} # number of vertices per element
{len(elements)} # number of elements
# Elements
"""
            )
            for element in elements:
                fobj.write(
                    " ".join(
                        _comsol_node_order(
                            [str(self.nodes[node]) for node in element], ecls
                        )
                    )
                    + " \n"
                )
            fobj.write(
                f"""
{len(elements)} # number of geometric entity indices
# Geometric entity indices
"""
            )
            fobj.write(
                " \n".join([str(self.elements[element]) for element in elements])
                + " \n"
            )

    def write_jcm(self, fobj):
        """Write a JCMsuite file (.jcm).

        Args:
            fobj (file object): file to write to
        """
        self.tag_jcm()
        elements_by_type = self.elements_by_type()
        binary = not isinstance(fobj, _io.TextIOBase)
        header = f"""/* <BLOBHead>
__BLOBTYPE__=Grid
__MODE__={'BINARY0' if binary else 'TEXT'}
__OPTIONS__=LHExchanged
__OWNER__=JCMwave
__TAG__={_random.randbytes(16).hex()}
__VERSION__=1.1.0
"""
        coordinate = self.info.get("jcm_header", {}).get("CoordinateSystem")
        if coordinate is not None:
            header += f"CoordinateSystem={coordinate}\n"
        if self._sdim == 2:
            element_types = (Triangle, Quadrilateral, Edge)
        elif self._sdim == 3:
            element_types = (
                Tetrahedron,
                Pyramid,
                Prism,
                Hexahedron,
                Triangle,
                Quadrilateral,
            )
        header += "\n".join(
            sorted(
                [
                    f"<I>ManifoldDim={max(self.dims)}",
                    f"<I>NPoints={len(self.nodes)}",
                    "<I>NRefinementSteps=0",
                    f"<I>SpaceDim={self._sdim}",
                ]
                + [
                    f"{ecls.jcm_info[self._sdim]}={len(elements_by_type.get(ecls, []))}"
                    for ecls in element_types
                ]
            )
        )
        cmin, cmax = _coord_minmax(self.nodes)
        header += f"""
<F>BBox_x_max={cmax[0]:.15e}
<F>BBox_x_min={cmin[0]:.15e}
<F>BBox_y_max={cmax[1]:.15e}
<F>BBox_y_min={cmin[1]:.15e}
<F>BBox_z_max={cmax[2]:.15e}
<F>BBox_z_min={cmin[2]:.15e}
<F>OriginX=0.000000000000000e+00
<F>OriginY=0.000000000000000e+00
<F>OriginZ=0.000000000000000e+00
<F>Twist=0.000000000000000e+00
*/
"""
        jcm_boundary = self.info.get("jcm_boundary", {})

        if binary:
            fobj.write(header.encode())
            fobj.write(b"0")
            for node, tag in self.nodes.items():
                fobj.write(_struct.pack(f"<i{self._sdim}d", tag, *node[: self._sdim]))
            for ecls in element_types:
                if ecls not in elements_by_type:
                    continue
                for element in elements_by_type[ecls]:
                    fobj.write(
                        _struct.pack(
                            f"<{element.n_nodes + 1}i",
                            *_jcm_node_order(
                                [self.nodes[node] for node in element], ecls
                            ),
                            self.elements[element],
                        )
                    )
                    if element.dim < self._sdim:
                        boundary = jcm_boundary.get(element, [1])
                        fobj.write(_struct.pack(f"{len(boundary)}i", *boundary))
            return

        fobj.write(header)
        fobj.write("# Points ###################\n# Points      (n, x, y, z)\n")
        for node, tag in self.nodes.items():
            fobj.write(
                f"# P\n{tag}\n"
                + "\n".join(f"{i:.15e}" for i in node[: self._sdim])
                + "\n"
            )
        boundary_header = True
        for ecls in element_types:
            if ecls not in elements_by_type:
                continue
            if boundary_header and ecls.dim < self._sdim:
                fobj.write("# Boundary patches ###################\n")
                boundary_header = False
            fobj.write(ecls.jcm_header[self._sdim])
            for element in elements_by_type[ecls]:
                fobj.write(
                    element.jcm_short_name
                    + "\n"
                    + "\n".join(
                        _jcm_node_order(
                            [str(self.nodes[node]) for node in element], ecls
                        )
                    )
                    + f"\n{self.elements[element]}\n"
                )
                if element.dim < self._sdim:
                    fobj.write(
                        "\n".join([str(i) for i in jcm_boundary.get(element, [1])])
                        + "\n"
                    )

    def write_gmsh(self, fobj):
        """Write a gmsh file (.msh).

        Args:
            fobj (file object): file to write to
        """
        self.tag_gmsh()
        fobj.write(
            """$MeshFormat
4.1 0 8
$EndMeshFormat
"""
        )
        physnames = self.info.get("physnames")
        if physnames is not None:

            fobj.write(f"$PhysicalNames\n{len(physnames)}\n")
            for (dim, tag), physname in physnames.items():
                fobj.write(f"{dim} {tag} {physname}\n")
            fobj.write("$EndPhysicalNames\n")
        domains_by_dim = dict(
            sorted(
                [((d.dim, tag), d) for d, tag in self.domains.items()],
                key=lambda item: item[0],
            )
        )

        fobj.write(
            f"""$Entities
{' '.join(str(sum(ddim == dim for ddim, _ in domains_by_dim)) for dim in range(4))}
"""
        )
        for (dim, tag), domain in domains_by_dim.items():
            if dim == 0:
                if len(domain) > 1:
                    raise ValueError("Entity of dimension 0 can only contain one Point")
                if len(domain) == 0:
                    continue
                x, y, z = next(iter(domain))[0].in_3d
                fobj.write(f"{tag} {x:.16g} {y:.16g} {z:.16g}")
            else:
                cmin, cmax = _coord_minmax(
                    {node for element in domain for node in element}
                )
                fobj.write(
                    f"{tag} {cmin[0]:.16g} {cmin[1]:.16g} {cmin[2]:.16g} "
                    f"{cmax[0]:.16g} {cmax[1]:.16g} {cmax[2]:.16g}"
                )

            phystags = self.info.get("phystags", {}).get(domain, [])
            fobj.write(f" {' '.join(map(str, [len(phystags)] + phystags))}")
            if dim == 0:
                fobj.write(" \n")
                continue
            bpoints = self.info.get("boundingpoints", {}).get(domain, [])
            fobj.write(f" {' '.join(map(str, [len(bpoints)] + bpoints))} \n")
        nodes = {node: tag for node, tag in self.nodes.items() if tag > 0}
        fobj.write(
            f"""$EndEntities
$Nodes
{len(self.domains)} {len(nodes)} {min(nodes.values())} {max(nodes.values())}
"""
        )
        for (dim, tag), domain in domains_by_dim.items():
            nodes_domain = {}
            for element in domain:
                for node in element:
                    nodetag = nodes.pop(node, None)
                    if nodetag is not None:
                        nodes_domain[node] = nodetag
            nodes_domain = dict(sorted(nodes_domain.items(), key=lambda item: item[1]))
            fobj.write(f"{dim} {tag} 0 {len(nodes_domain)}\n")
            for nodetag in nodes_domain.values():
                fobj.write(f"{nodetag}\n")
            for node in nodes_domain:
                x, y, z = node.in_3d
                fobj.write(f"{x:.16g} {y:.16g} {z:.16g}\n")
        fobj.write("$EndNodes\n$Elements\n")
        nblocks = sum(
            len({type(e) for e in domain}) for domain in self.domains if domain.dim > 0
        )
        nelements = sum(len(domain) for domain in self.domains if domain.dim > 0)

        fobj.write(
            f"{nblocks} {nelements} "
            f"{min(self.elements.values())} {max(self.elements.values())}\n"
        )
        elements_by_domain_and_type = {}

        for (dim, _), domain in domains_by_dim.items():
            val = np.array([d for d in domain], dtype=object)
            vallist = [d for d in domain]
            ds = []
            gmv = np.array([d.gmsh_type for d in domain])
            gsort = np.sort(gmv)
            valsort = val[np.argsort(gmv)]
            uni, ind = np.unique(gsort, return_index=True)
            if dim == 0:
                continue
            elements_by_domain_and_type[domain] = {}
            groups = np.split(valsort, ind[1:])
            lsgroup = []
            ind = np.append(ind, len(val))
            for i in range(len(ind) - 1):
                lsgroup.append(vallist[ind[i] : ind[i + 1]])
            for i, u in enumerate(uni):
                elements_by_domain_and_type[domain][u] = lsgroup[i]
        gmsh_types = sorted(
            [e.gmsh_type for e in Element.__subclasses__() if e.gmsh_type is not None]
        )
        for domain, elements_by_type in elements_by_domain_and_type.items():
            domaintag = self.domains[domain]
            for gmsh_type in gmsh_types:
                elements = elements_by_type.get(gmsh_type, None)
                if elements is None:
                    continue

                fobj.write(f"{domain.dim} {domaintag} {gmsh_type} {len(elements)}\n")
                for element in elements:
                    fobj.write(
                        f"{self.elements[element]} "
                        f"{' '.join(str(self.nodes[node]) for node in element)} \n"
                    )
        fobj.write("$EndElements\n")

    def gmsh(self, model):
        """Create a gmsh model.

        Add the stored mesh to the model given. The model has to be initialized before
        handing it to this method. Afterwards, it can be attempted to generate surfaces
        with `model.mesh.createTopology()`, if they are not included in the mesh.

        Args:
            model (gmsh.Model): gmsh model

        Returns:
            gmsh.Model
        """
        self.tag_gmsh()
        for domain, tag in self.domains.items():
            model.addDiscreteEntity(domain.dim, tag)

            nodes_coords = []
            nodes_tags = []
            for element in domain:
                for node in element:
                    if self.nodes[node] not in nodes_tags:
                        nodes_coords.extend(node.in_3d)
                        nodes_tags.append(self.nodes[node])
            model.mesh.addNodes(domain.dim, tag, nodes_tags, nodes_coords)

            element_tags = {}
            element_nodes = {}
            for element in domain:
                element_tags.setdefault(element.gmsh_type, []).append(
                    self.elements[element]
                )
                element_nodes.setdefault(element.gmsh_type, []).extend(
                    [self.nodes[node] for node in element]
                )
            model.mesh.addElements(
                domain.dim,
                tag,
                list(element_tags.keys()),
                list(element_tags.values()),
                list(element_nodes.values()),
            )
        return model


read = Mesh.read


# @jit(nopython=True)
def loopjit(
    ed,
):
    for element in domain:

        ed.setdefault(element.gmsh_type, []).append(element)


def _read_header_jcm(fobj):
    """Read the header of a JCMsuite file.

    It is assumed that the header starts at the first line of the given file descriptor.
    The start of the header is defined by '/* <BLOBHead>' and its end by the first
    occurence of '*/'. It reads the coordinate system info from the header and returns
    it in a dict if found.

    Args:
        fobj (file object): handle to the beginning of a jcm file

    Returns:
        dict
    """
    line = next(fobj)
    binary = False
    if isinstance(line, bytes):
        line = line.decode()
        binary = True
    if line != "/* <BLOBHead>\n":
        raise ValueError("missing header in jcm file")
    header = {}
    for line in fobj:
        line = line.decode() if binary else line
        if line == "*/\n":
            break
        key, *val = line.split("=")
        val = "=".join(val)
        header[key] = val[:-1]
    return header


def _read_jcm(handle, mode, binary=False):
    """Read numbers from JCM file.

    The next numbers for the file are read and convert according to mode. The file is
    read either as a text or binary file. If the file is in text mode, each line must
    contain exactly one number.

    Args:
        handle (file object): file in text or binary mode
        mode (str): mode definition similar to `struct.unpack`, `idd` means one integer
            and two doubles and can also be written `i2d`
        binary (bool, default False): read file as text or binary file

    Returns:
        tuple
    """
    ndigits = 0
    if binary:
        vals = {"i": 4, "d": 8}
        nbytes = 0
        for char in mode:
            if char.isdigit():
                ndigits = ndigits * 10 + int(char)
                continue
            nbytes += max(ndigits, 1) * vals[char]
            ndigits = 0
        return _struct.unpack("<" + mode, handle.read(nbytes))
    vals = {"i": int, "d": float}
    res = []
    for char in mode:
        if char.isdigit():
            ndigits = ndigits * 10 + int(char)
            continue
        for _ in range(max(ndigits, 1)):
            res.append(vals[char](next(handle)))
        ndigits = 0
    return tuple(res)


def _iter_wo_comments(fobj, indicator="#"):
    """Remove comments.

    This is a crude method to remove comments, that is unaware of strings, escape
    characters or anything else.

    Args:
        fobj (file object): file
        indicator (str, default '#'): indicator for the start of a comment

    Yields:
        str
    """
    for line in fobj:
        if isinstance(line, bytes):
            line = line.decode()
        line = line.split(indicator)[0].strip()
        if line:
            yield line


def _comsol_node_order(nodes, element_cls):
    """Change the node order to/from COMSOL style.

    It is exploited that COMSOL just swaps the postions of two nodes and not more
    complicated reorderings.

    Args:
        nodes (list): nodes of the element
        element_cls (object): type of element from which the nodes are

    Returns:
        list
    """
    if element_cls == Quadrilateral:
        return nodes[0], nodes[1], nodes[3], nodes[2]
    if element_cls == Pyramid:
        return nodes[0], nodes[1], nodes[3], nodes[2], nodes[4]
    if element_cls == Hexahedron:
        return (
            nodes[0],
            nodes[1],
            nodes[3],
            nodes[2],
            nodes[4],
            nodes[5],
            nodes[7],
            nodes[6],
        )
    return nodes


def _jcm_node_order(nodes, element_cls):
    """Change the node order to/from JCMsuite style.

    It is exploited that JCMsuite just swaps the postions of two nodes and not more
    complicated reorderings.

    Args:
        nodes (list): nodes of the element
        element_cls (object): type of element from which the nodes are

    Returns:
        List
    """
    if element_cls == Hexahedron:
        return (
            nodes[0],
            nodes[1],
            nodes[3],
            nodes[2],
            nodes[4],
            nodes[5],
            nodes[7],
            nodes[6],
        )
    return nodes


def _gmsh_iter(fobj):
    """Iterater through gmsh file.

    Ignores all sections except for a list of known ones.

    Args:
        fobj (file object): file

    Yields:
        str
    """
    modes = [
        "$MeshFormat",
        "$PhysicalNames",
        "$Entities",
        "$Nodes",
        "$Elements",
    ]
    itr = _iter_wo_comments(fobj, "//")
    mode = None
    for line in itr:
        if mode is None and line.startswith("$"):
            mode = line
        if mode in modes:
            yield line
        if mode is not None and line == "$End" + mode[1:]:
            mode = None


def _coord_minmax(nodes):
    """Calculate the maximal and minimal coordinate value for each direction.

    If the list is empty return zeros.

    Args:
        nodes (Iterable): list of nodes

    Returns:
        tuple
    """
    if not nodes:
        return [0, 0, 0], [0, 0, 0]
    cmin = [float("inf")] * 3
    cmax = [float("-inf")] * 3
    for node in nodes:
        for i in range(3):
            cmin[i] = min((node.in_3d[i], cmin[i]))
            cmax[i] = max((node.in_3d[i], cmax[i]))
    return cmin, cmax


def _cli():
    """Entry point for the CLI."""
    parser = _argparse.ArgumentParser(
        prog="meshparser",
        description="convert between comsol, jcmsuite, and gmsh meshes",
    )
    parser.add_argument("input", help="with extension jcm, mphtxt, or msh")
    parser.add_argument("output", help="with extension jcm, mphtxt, or msh")
    parser.add_argument(
        "-m",
        "--mode",
        help="writing mode (binary only possible for .jcm files)",
        choices=["a", "w", "x", "ab", "wb", "xb"],
        default="x",
    )
    parser.add_argument(
        "-s",
        "--skip_1d",
        help="skip writing 1D elements in 2D COMSOL meshes",
        action="store_true",
    )
    # JUAT COMMENTED THIS SORRY
    # parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()
    mesh = read(args.input)
    mesh.write(args.output, mode=args.mode, skip_1d=args.skip_1d)
    return 0


if __name__ == "__main__":
    _cli()
