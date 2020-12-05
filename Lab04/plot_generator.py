import matplotlib.pyplot as plt
import numpy as np

def plot_map(x, y, z, omega, title='', zlabel=''):

    plt.xlabel("x")
    plt.ylabel("y")

    plt.title(title + ' relaxation - ' + zlabel + ' - omega=' + str(omega), fontweight='bold')
    figure = plt.gcf()
    figure.set_size_inches(12, 8)
    z_min, z_max = np.amin(z), np.amax(z)
    plt.pcolor(x, y, z, cmap='seismic', vmin=z_min, vmax=z_max)
    plt.colorbar()
    
    plt.savefig("{}_omega-{}_{}_relaxation.png".format(zlabel, str(omega), title))
    # plt.show()
    plt.close()

def plot(x, y, omega, title='', xlabel='', ylabel=''):

    for i, x_current in enumerate(x):
        label = 'omega=' + str(omega[i]) + '; it=' + str(x_current[-1])
        plt.plot(x_current, y[i], label=label)

    plt.xscale('log')
    plt.xlim(left=1)

    plt.title(title + ' relaxation - ' + ylabel, fontweight='bold')
    plt.xlabel(xlabel, fontweight='bold')
    plt.ylabel(ylabel, fontweight='bold')

    plt.grid(True)
    plt.legend(loc="upper right")
    plt.savefig("{}_{}_relaxation.png".format(ylabel, title))
    # plt.show()
    plt.close()