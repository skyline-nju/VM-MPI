import numpy as np
import matplotlib.pyplot as plt
import struct
import os


def read_rhox(fin):
    """ read the density profile averaged over y direction, rho_y(x, t)."""

    basename = os.path.basename(fin)
    str_list = basename.rstrip(".bin").lstrip("rhox_").split(".")
    eta = float(str_list[0]) / 1000
    eps = float(str_list[1]) / 1000
    Lx = int(str_list[2])
    Ly = int(str_list[3])
    frame_size = 4 * Lx
    with open(fin, "rb") as f:
        f.seek(0, 2)
        file_size = f.tell()
        f.seek(0, 0)
        while f.tell() < file_size:
            buf = f.read(frame_size)
            rhox = np.array(struct.unpack("%df" % Lx, buf))
            yield eta, eps, Lx, Ly, rhox


def read_filed(fin, lbox=4):
    """ read coarse-grained fields for rho, vx, vy. """
    basename = os.path.basename(fin)
    str_list = basename.rstrip(".bin").split("_")
    L = int(str_list[2])
    n = L // lbox
    frame_size = n * n * 12
    with open(fin, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()
        f.seek(0)
        while f.tell() < filesize:
            buf = f.read(frame_size)
            data = np.array(struct.unpack("%df" % (n * n * 3), buf))
            rho, vx, vy = data.reshape(3, n, n)
            yield rho, vx, vy


def read_snap(fin):
    """ read snapshot saved in the form {x_1, y_1, theta_1, ..., x_i, y_i,
        theta_i, ...}
    """
    basename = os.path.basename(fin)
    str_list = basename.rstrip("bin").lstrip("s").split(".")
    # eta = float(str_list[0]) / 1000
    # eps = float(str_list[1]) / 1000
    Lx = int(str_list[2])
    Ly = int(str_list[3])
    with open(fin, "rb") as f:
        buf = f.read()
        data = struct.unpack("%df" % (Lx * Ly * 3), buf)
        x, y, theta = np.array(data).reshape(Lx * Ly, 3).T
    return x, y, theta
    

if __name__ == "__main__":
    rhox_file = r"data\rhox_350.20.300.100.1.bin"
    field_file = r"data\RT_ave_512_0.180_0.000_200_1.bin"
    frames = read_filed(field_file)
    for (rho, vx, vy) in frames:
        rho_m = np.mean(rho)
        print(rho_m)
        # vxm = np.mean(vx) / rho_m
        # vym = np.mean(vy) / rho_m
        # print(np.sqrt(vxm ** 2 + vym ** 2))

    
    # snap_file = r"VM_MPI_v2\data\s350.0.300.100.1.0004.bin"
    # x, y, theta = read_snap(snap_file)
    # plt.plot(x, y, ".")
    # plt.show()
