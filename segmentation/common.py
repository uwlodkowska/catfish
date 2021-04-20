#!/usr/bin/env python
# coding: utf-8

class Approximation(): 
    def __init__(self, reg_props, approx_patch, intensity,  idno):
        self.reg_props = reg_props
        self.approx_patch = approx_patch
        self.intensity = intensity
        self.idno = idno
        


def save_approximations_overlay(img, found, fig_path):
    ax = plt.subplot(aspect='equal')
    plt.imshow(img)
    for instance in found:
        ax.add_artist(instance.approx_patch)
    figure = plt.gcf() # get current figure
    figure.set_size_inches(img.shape[1]/100, img.shape[0]/100)
    plt.savefig(fig_path, dpi=100)
    plt.clf()
