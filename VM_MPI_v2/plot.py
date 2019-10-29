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
    rhox_file = r"VM_MPI_v2\data\rhox_350.20.300.100.1.bin"
    frames = read_rhox(rhox_file)

    for eta, eps, Lx, Ly, rhox in frames:
        plt.plot(rhox)
        plt.show()
        plt.close()
    
    # snap_file = r"VM_MPI_v2\data\s350.0.300.100.1.0004.bin"
    # x, y, theta = read_snap(snap_file)
    # plt.plot(x, y, ".")
    # plt.show()
