#!/usr/bin/env python3
#! -*- coding: utf-8 -*-

"""
Get transformed fractions and hardness for different cooling rates 
"""
import argparse
import matplotlib.pyplot as plt
from transformation_models import Alloy, TransformationDiagrams

if __name__ == '__main__':
    # Defines alloy (grain size gs and composition)
    alloy = Alloy(gs=7, C=0.37, Mn=0.77, Si=0.15, Ni=0.04, Cr=0.98, Mo=0.21)

    # Initializes diagrams object
    diagrams = TransformationDiagrams(alloy)

    Tini = 900.  # initial temperature
    Tfin = 25.  # final temperature
    cooling_rates = [1000, 300, 100, 30, 10, 3, 1, 3e-1,
                     1e-1, 3e-2, 1e-2, 3e-3, 1e-3]  # cooling rates, duh

    f_ferr = []  # list with fractions of ferrite at Tfin
    f_pear = []  # list with fractions of pearlite at Tfin
    f_bain = []  # list with fractions of bainite at Tfin
    f_mart = []  # list with fractions of martensite at Tfin
    Hv = []  # list with hadness values at Tfin

    for phi in cooling_rates:
        print('Calculating phase fractions for phi={:g} oC/s...'.format(phi))

        total_time = (Tini - Tfin)/phi  # Total heat treatment time
        # Gets phase fractions and hardness
        f = diagrams.get_transformed_fraction([0, total_time], [Tini, Tfin])

        f_fin = f.iloc[-1]  # f at Tfin

        # Appends phase fractions and hardness at Tfin
        f_ferr.append(f_fin['pearlite'])
        f_pear.append(f_fin['ferrite'])
        f_bain.append(f_fin['bainite'])
        f_mart.append(f_fin['martensite'])
        Hv.append(f_fin['Hv'])

    # Initializes plot window
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(wspace=.2)

    # Plot
    ax1.plot(cooling_rates, f_ferr, label='Pearlite')
    ax1.plot(cooling_rates, f_pear, label='Ferrite')
    ax1.plot(cooling_rates, f_bain, label='Bainite')
    ax1.plot(cooling_rates, f_mart, label='Martensite')
    ax1.set_xlabel(u'Cooling rate (°C/s)')
    ax1.set_ylabel(u'Phase fraction')
    ax1.set_xscale('log')
    ax1.set_title(u'Phase fractions at {:g} °C'.format(Tfin))
    ax1.legend()

    ax2.plot(cooling_rates, Hv)
    ax2.set_xlabel(u'Cooling rate (°C/s)')
    ax2.set_ylabel(u'Vickers Hardness')
    ax2.set_title(u'Hardness for phase fractions at {:g} °C'.format(Tfin))
    ax2.set_xscale('log')

    fig.suptitle(alloy.format_composition())

    plt.show()
