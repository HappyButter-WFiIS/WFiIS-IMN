import matplotlib.pyplot as plt

def plot(x, y, title='', xlabel='', ylabel=''):

    plt.plot(x[0], y[0], label='tol_1=10^-2')
    plt.plot(x[1], y[1], label='tol_2=10^-5')

    plt.title(ylabel+ ' ' + title + ' method', fontweight='bold')
    plt.xlabel(xlabel, fontweight='bold')
    plt.ylabel(ylabel, fontweight='bold')

    plt.grid(True)
    plt.legend(loc="upper right")
    plt.savefig("{}_{}.png".format(ylabel,title))
    plt.close()