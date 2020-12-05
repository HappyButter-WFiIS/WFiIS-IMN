import matplotlib.pyplot as plt
import numpy as np

def plot(t, u, title='', xlabel='', ylabel=''):
    N = 500

    plt.plot(t, u, label='u(t)')
    z = N - u
    plt.plot(t, z, label='z(t)')

    plt.title(title + ' method', fontweight='bold')
    plt.xlabel(xlabel, fontweight='bold')
    plt.ylabel(ylabel, fontweight='bold')

    plt.xticks(np.arange(0,101,10))
    plt.yticks(np.arange(0,501,50))
    plt.grid(True)

    plt.savefig("{}.png".format(title))
    plt.close()
