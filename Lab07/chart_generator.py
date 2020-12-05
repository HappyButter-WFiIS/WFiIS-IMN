import matplotlib.pyplot as plt
import numpy as np


def plot_map(x, y, z, Q, title=''):

    plt.xlabel("x")
    plt.ylabel("y")

    plt.title('Q=' + str(Q) + ', ' + title, fontweight='bold')

    figure = plt.gcf()
    figure.set_size_inches(12, 8)
    
    z_min, z_max = np.amin(z), np.amax(z)
    plt.pcolor(x, y, z, cmap='afmhot', vmin=z_min, vmax=z_max)
    plt.colorbar()
    
    plt.savefig("Q_{}-{}.png".format(str(Q), title))
    # plt.show()
    plt.close()


def plot_contour_map(x, y, z, Q, title=''):
    plt.xlabel("x")
    plt.ylabel("y")

    plt.title('Q=' + str(Q) + ', ' + title, fontweight='bold')

    figure = plt.gcf()
    figure.set_size_inches(12, 8)

    custom_cmap = 'afmhot'
    if(Q==4000):
        custom_cmap = 'PiYG'

    plt.contour(x, y, z, cmap=custom_cmap,levels=40)
    plt.colorbar()

    plt.savefig("Q_{}-{}.png".format(str(Q), title))
    # plt.show()
    plt.close()