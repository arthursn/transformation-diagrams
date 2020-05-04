#!/usr/bin/env python3
#! -*- coding: utf-8 -*-

"""
Plot TTT (or CCT) diagram and transformed fraction
"""
import argparse
import matplotlib.pyplot as plt
from transformation_models import Alloy, TransformationDiagrams

if __name__ == '__main__':
    # Defines alloy (grain size gs and composition)
    alloy = Alloy(gs=7, C=0.37, Mn=0.77, Si=0.15, Ni=0.04, Cr=0.98, Mo=0.21)

    # Initializes diagrams object
    diagrams = TransformationDiagrams(alloy)

    """    
    Example 1: Plot TTT and CCT diagrams
    """
    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig1.subplots_adjust(wspace=.2)

    diagrams.TTT(ax=ax1)
    ax1.set_xlim(1e-2, 1e8)
    ax1.set_ylim(300, 1000)

    diagrams.CCT(ax=ax2)
    ax2.set_xlim(1e-2, 1e8)
    ax2.set_ylim(300, 1000)

    fig1.suptitle(ax1.get_title())
    ax1.set_title('')
    ax2.set_title('')

    """
    Example 2: Get phase fraction information and save data in file
    """
    # Chooses a thermal cycle (continuous cooling from 1000 to 0 oC in
    # 2000 s), draws cooling curve in the CCT and plots phase fraction
    t, T = [0, 2000], [1000, 0]

    # get_transformed_fraction returns a pandas DataFrame
    # n is number of points
    data = diagrams.get_transformed_fraction(t, T, n=2001)

    # Displays data
    print(data)

    # Save as csv
    data.to_csv('data.csv', index=False)

    # Save as excel
    data.to_excel('data.xls', index=False)

    """
    Example 3: Plot CCT diagram and transformed fraction
    """
    fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 6))
    fig2.subplots_adjust(wspace=.2)

    # Plot CCT diagram
    diagrams.CCT(Tini=1000, ax=ax3)

    # Same thermal cycle as from example 2
    diagrams.draw_thermal_cycle(ax3, t, T)
    diagrams.plot_phase_fraction(t, T, xaxis='T', ax=ax4)

    fig2.suptitle(ax3.get_title())
    ax3.set_title('')
    ax4.set_title('')

    plt.show()
