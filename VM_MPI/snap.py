import struct
import matplotlib.pyplot as plt
import numpy as np


def read(file):
    with open(file, "rb") as f:
        buff = f.read()
        print(len(buff))
        data = struct.unpack("%dB" % (108000 * 10), buff)
        data = np.array(data).reshape(10, 600, 180)
        return data


if __name__ == "__main__":
    file = "c_0.35_0.02_180_600_108000_1_100.bin"
    data = read(file)
    for frame in data:
        plt.imshow(frame.T, origin="lower")
        plt.show()
        plt.close()
