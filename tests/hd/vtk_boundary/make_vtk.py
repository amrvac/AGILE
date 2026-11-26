"""Write the two boundary-data files of this test (legacy ASCII VTK).

mixed.vtk      face min1, density only, three times (time-varying)
primitive.vtk  face min3, all five primitive variables, one time (static)

Axes are the face's two tangential coordinates, then time; the first varies
fastest. Run from this directory: python make_vtk.py
"""
import math

N = 9                      # samples per tangential direction over [0, 1]
H = 1.0 / (N - 1)


def write(name, title, nt, dt, fields):
    with open(name, "w") as f:
        f.write("# vtk DataFile Version 2.0\n%s\nASCII\n" % title)
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write("DIMENSIONS %d %d %d\n" % (N, N, nt))
        f.write("SPACING %r %r %r\n" % (H, H, dt))
        f.write("ORIGIN 0 0 0\n")
        f.write("POINT_DATA %d\n" % (N * N * nt))
        for vname, fun in fields:
            f.write("SCALARS %s double 1\nLOOKUP_TABLE default\n" % vname)
            for k in range(nt):
                for j in range(N):
                    f.write(" ".join("%.17g" % fun(i * H, j * H, k * dt)
                                     for i in range(N)) + "\n")


def bump(a, b):
    return math.sin(math.pi * a) * math.sin(math.pi * b)


write("mixed.vtk", "hd density on face min1", 3, 0.1,
      [("rho", lambda y, z, t: 1.0 + (0.5 + 2.0 * t) * bump(y, z))])

write("primitive.vtk", "hd primitive state on face min3", 1, 1.0,
      [("rho", lambda x, y, t: 1.0 + 0.5 * bump(x, y)),
       ("v1", lambda x, y, t: 0.1 * (x - 0.5)),
       ("v2", lambda x, y, t: -0.1 * (y - 0.5)),
       ("v3", lambda x, y, t: 0.5 + 0.2 * bump(x, y)),
       ("p", lambda x, y, t: 1.0 + 0.3 * bump(x, y))])
