import matplotlib.pyplot as plt
import numpy as np

def plot_map(x, y, z, k, title='', zlabel=''):

    plt.xlabel("x")
    plt.ylabel("y")

    plt.title(title + ' wielosiatkowa - ' + zlabel + ' - k=' + str(k), fontweight='bold')
    figure = plt.gcf()
    figure.set_size_inches(12, 8)
    z_min, z_max = np.amin(z), np.amax(z)

    plt.pcolor(x, y, z, cmap='seismic', vmin=z_min, vmax=z_max)
    plt.colorbar()
    
    plt.savefig("{}_k={}_{}_wielosiatkowa.png".format(zlabel, str(k), title))
    # plt.show()
    plt.close()

def plot(k_arr, s_arr, title='', xlabel='', ylabel=''):

    for i, s in enumerate(s_arr):
        it_max = len(s)
        it_array = np.linspace(1, it_max, it_max, endpoint=True, dtype=np.int16)

        label = "K = " + str(k_arr[i])

        plt.plot(it_array, s, label=label)

    plt.title(title + ' wielosiatkowa - ' + ylabel, fontweight='bold')
    plt.xlabel(xlabel, fontweight='bold')
    plt.ylabel(ylabel, fontweight='bold')

    plt.grid(True)
    plt.legend(loc="upper right")
    plt.savefig("{}_{}_wielosiatkowa.png".format(ylabel, title))
    # plt.show()
    plt.close()