import matplotlib.pyplot as plt
import numpy as np

def generate_plot_area(title1='', title2='', xlabel1='', xlabel2='', ylabel1='', ylabel2=''):

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(20, 8)
    ax1.set_title(title1, fontweight='bold')
    ax1.set_xlabel(xlabel1, fontweight='bold')
    ax1.set_ylabel(ylabel1, fontweight='bold')

    ax2.set_title(title2, fontweight='bold')
    ax2.set_xlabel(xlabel2, fontweight='bold')
    ax2.set_ylabel(ylabel2, fontweight='bold')

    return ax1, ax2
