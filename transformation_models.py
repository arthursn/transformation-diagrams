#!/usr/bin/env python3
#! -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from abc import abstractmethod
from scipy import integrate
from scipy.interpolate import splrep, splev, interp1d
from scipy.optimize import root

R = 8.314459
K = 273.15


def parse_comp(**comp):
    """
    Parse a composition dictionary and returns default values if
    values not set
    """
    C = comp.pop('C', 0)
    Mn = comp.pop('Mn', 0)
    Si = comp.pop('Si', 0)
    Ni = comp.pop('Ni', 0)
    Cr = comp.pop('Cr', 0)
    Mo = comp.pop('Mo', 0)
    Co = comp.pop('Co', 0)
    return C, Mn, Si, Ni, Cr, Mo, Co


def FahrenheitToCelsius(TF):
    return (TF - 32.)*5./9.


def FC(**comp):
    """
    Function that expresses the effects of the alloying elements on
    on the kinetics of ferrite transformation
    """
    C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**comp)
    return np.exp((1.0 + 6.31*C + 1.78*Mn + 0.31*Si + 1.12*Ni + 2.7*Cr + 4.06*Mo))


def PC(**comp):
    """
    Function that expresses the effects of the alloying elements on
    on the kinetics of pearlite transformation
    """
    C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**comp)
    return np.exp(-4.25 + 4.12*C + 4.36*Mn + 0.44*Si + 1.71*Ni + 3.33*Cr + 5.19*np.sqrt(Mo))


def BC(**comp):
    """
    Function that expresses the effects of the alloying elements on
    on the kinetics of bainite transformation
    """
    C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**comp)
    return np.exp(-10.23 + 10.18*C + 0.85*Mn + 0.55*Ni + 0.9*Cr + 0.36*Mo)


def Ae1_Grange(**comp):
    """
    Grange's equation for Ae1
    """
    C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**comp)
    return FahrenheitToCelsius(1333 - 25*Mn + 40*Si - 26*Ni - 42*Cr)


def Ae3_Grange(**comp):
    """
    Grange's equation for Ae3
    """
    C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**comp)
    return FahrenheitToCelsius(1570 - 323*C - 25*Mn + 80*Si - 32*Ni - 3*Cr)


def Ae1_Andrews(**comp):
    """
    Andrews' equation for Ae1
    """
    C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**comp)
    return 712 - 17.8*Mn + 20.1*Si - 19.1*Ni + 11.9*Cr + 9.8*Mo


def Ae3_Andrews(**comp):
    """
    Andrews' equation for Ae3
    """
    C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**comp)
    return 871 - 254.4*np.sqrt(C) + 51.7*Si - 14.2*Ni


def Bs_Li(**comp):
    """
    Bs calculation from Li's work
    """
    C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**comp)
    return 637 - 58*C - 35*Mn - 15*Ni - 34*Cr - 41*Mo


def Bs_VanBohemen(**comp):
    """
    [1] S.M.C. van Bohemen, Mater. Sci. Technol. 28 (2012) 487–495.
    """
    C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**comp)
    return 839 - (86*Mn + 23*Si + 67*Cr + 33*Ni + 75*Mo) - 270*(1 - np.exp(-1.33*C))


def Ms_Andrews(**comp):
    """
    Andrews' equation for Ms
    """
    C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**comp)
    return 539 - 423*C - 30.4*Mn - 17.7*Ni - 12.1*Cr - 7.5*Mo + 10*Co - 7.5*Si


def Ms_VanBohemen(**comp):
    """
    [1] S.M.C. van Bohemen, Mater. Sci. Technol. 28 (2012) 487–495.
    """
    C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**comp)
    return 565 - (31*Mn + 13*Si + 10*Cr + 18*Ni + 12*Mo) - 600*(1-np.exp(-0.96*C))


class Alloy:
    """
    Alloy properties (composition in wt.% and prior austenite grain size)
    """

    def __init__(self, gs, **w):
        # Grain size
        self.gs = gs

        # Alloy composition
        self.w = w
        # Main elements
        self.C, self.Mn, self.Si, self.Ni, self.Cr, self.Mo, self.Co = parse_comp(**w)

        self.FC = FC(**w)
        self.PC = PC(**w)
        self.BC = BC(**w)
        # self.Ae3 = Ae3_Andrews(**w)
        # self.Ae1 = Ae1_Andrews(**w)
        self.Ae3 = Ae3_Grange(**w)
        self.Ae1 = Ae1_Grange(**w)
        # self.Ae3 = 970
        # self.Ae1 = 590
        self.Bs = Bs_Li(**w)
        # self.Ms = Ms_Andrews(**w)
        # self.Bs = Bs_VanBohemen(**w)
        self.Ms = Ms_VanBohemen(**w)

    def format_composition(self, vmin=0):
        fmt = []
        for k, v in self.w.items():
            if v > vmin:
                fmt.append('{:g}{:}'.format(v, k))
        fmt.insert(0, 'Fe')  # assumes it is steel
        return '-'.join(fmt) + ' (wt.%)'


class SigmoidalFunction(object):
    """
    Abstract class for S(X) and I(X) functions. Once initialized,
    calculates values of the function for a given [xmin, xmax]
    interval and then creates a spline interpolator. The returned
    values are calculated by the interpolator. This method has the
    advantage of being able to process x as an array (or any other
    kind of iterator)
    """
    tck = None  # tck spline knots, coefficients and degree
    tck_inv = None  # spline parameters of the inverse function

    def __new__(cls, x):
        """
        __new__ behaviour is modified to return the interpolated
        value of the function
        """
        if cls is SigmoidalFunction:
            raise TypeError("Can't instantiate abstract class SigmoidalFunction")

        # Check for compulsory subclass attributes
        for var in ['xmin', 'xmax', 'ymin', 'ymax', 'n']:
            if not hasattr(cls, var):
                raise NotImplementedError('Class {} lacks required `{}` class attribute'.format(cls, var))

        # This is were S(X) or I(X) is returned
        return cls.val(x)

    @staticmethod
    def f(x):
        """
        Function to be integrated
        """
        pass

    @classmethod
    def val(cls, x):
        """
        Evaluates SigmoidalFunction(x)
        """
        if hasattr(x, '__iter__') and not isinstance(x, str):
            x = np.array(x)
            xmin = x[x > 0].min()
        else:
            xmin = x
        xmin = min(cls.xmin, xmin)

        # init spline if not initialized yet or if xmin is lower than lower bound
        if xmin < cls.xmin or cls.tck is None:
            cls.xmin = xmin
            cls.init_spline()

        return splev(x, cls.tck)

    @classmethod
    def inv(cls, y):
        """
        Evaluates inverse function SigmoidalFunction^-1(y)
        """
        if hasattr(y, '__iter__') and not isinstance(y, str):
            y = np.array(y)
            ymin, ymax = y.min(), y.max()
        else:
            ymin, ymax = y, y

        if cls.tck_inv is None:
            cls.init_spline()

        if ymin < cls.ymin or ymax > cls.ymax:
            print('Be careful! y value out of bounds [{:g}:{:g}]. '
                  'Returned value is an extrapolation'.format(cls.ymin, cls.ymax))

        return splev(y, cls.tck_inv)

    @classmethod
    def init_spline(cls):
        """
        Initializes spline
        """
        X = np.linspace(cls.xmin, cls.xmax, cls.n)
        Y = np.array([integrate.quad(cls.f, 0, x)[0] for x in X])
        cls.ymin = Y.min()
        cls.ymax = Y.max()
        cls.tck = splrep(X, Y)
        cls.tck_inv = splrep(Y, X)


class S(SigmoidalFunction):
    n = 999
    xmin = 0.001
    xmax = 0.999
    ymin = 0.02638507
    ymax = 2.02537893

    """
    S(X) function calculated using numerical integration and
    spline interpolation
    """
    @staticmethod
    def f(x):
        return 1./(x**(0.4*(1. - x))*(1. - x)**(0.4*x))


class I(SigmoidalFunction):
    """
    I(X) function calculated using numerical integration and
    spline interpolation
    """
    n = 999
    xmin = 0.001
    xmax = 0.999
    ymin = 0.29961765
    ymax = 4.05928646

    @staticmethod
    def f(x):
        return 1./(x**(2.*(1. - x)/3.)*(1. - x)**(2.*x/3.))


class PhaseTransformation(object):
    """
    Abstract class for calculating kinetics of diffusional phase
    transformations
    """

    def __init__(self, alloy):
        self.alloy = alloy
        self.initialize()
        # Check for compulsory object attributes
        for var in ['comp_factor', 'Ts', 'Tf']:
            if not hasattr(self, var):
                raise NotImplementedError('Object {} lacks required `{}` attribute'.format(self, var))

    @classmethod
    def __init_subclass__(cls):
        # Check for compulsory subclass attributes
        for var in ['Q', 'n1', 'n2']:
            if not hasattr(cls, var):
                raise NotImplementedError('Class {} lacks required `{}` class attribute'.format(cls, var))

    @abstractmethod
    def initialize(self):
        pass

    def get_transformation_factor(self, T):
        return self.comp_factor/(2**(self.n1*self.alloy.gs)*(self.Ts - T)**self.n2*np.exp(-self.Q/(R*(T + K))))

    def get_transformation_time(self, T, f):
        return S(f)*self.get_transformation_factor(T)

    def get_transformation_temperature(self, Tmax, Tmin, cooling_rate, f, dT=1.0):
        """
        cooling_rate : iterable
        """
        dt = dT/np.array(cooling_rate)
        nt = len(dt) if hasattr(dt, '__len__') else 1
        T = np.arange(Tmax, Tmin, -dT)
        nucleation_time = np.full((nt, len(T)), 0, dtype=float)

        filtr = T < self.Ts
        nucleation_time[:, filtr] = np.outer(dt, 1./self.get_transformation_factor(T[filtr]))
        nucleation_time = nucleation_time.cumsum(axis=1)

        Tt = np.full(nt, np.nan, dtype=float)  # Transformation temperature for a given fraction f

        # Check for indices of nucleation_time larger than threshold S(f)
        # First occurrence is the transformation temperature
        Sf = S(f)
        for i, n_time in enumerate(nucleation_time):
            idx, = np.where(n_time >= Sf)
            if len(idx) > 0:
                Tt[i] = T[idx[0]]

        return float(Tt) if nt == 1 else Tt

    def get_transformed_fraction(self, t, T, n=1000):
        """
        t, T : iterables
        """
        if len(t) > 3:
            # Fits T(t) by spline
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            # Uses linear interpolator
            t2T = interp1d(t, T)

        # To ensure convergence of the algorithm, the T(t) thermal cycle is
        # adjusted by a spline and the nucleation time is calculated by
        # increments dt = (max(t) - min(t))/n
        dt = (max(t) - min(t))/(n - 1)
        t = np.linspace(min(t), max(t), n)
        T = t2T(t)
        nucleation_time = np.full(t.shape, 0, dtype=float)
        f = np.full(T.shape, 0, dtype=float)

        # Calculates nucleation time only for T lower than transformation
        # start temperature and higher than Tf
        filtr = (T < self.Ts) & (T > self.Tf)
        if np.any(filtr):
            nucleation_time[filtr] = dt/self.get_transformation_factor(T[filtr])
            nucleation_time = nucleation_time.cumsum()
            if T[0] < self.Ts:
                # This is the factor corresponding to the transformed fraction at t[0]
                nucleation_time += min(t)/self.get_transformation_factor(T[0])

            # New filter: calculates f only for nucleation_time inside the bounds
            # of S.inv(y)
            filtr = (nucleation_time >= S.ymin) & (nucleation_time <= S.ymax)
            if np.any(filtr):
                f[filtr] = S.inv(nucleation_time[filtr])
                f[nucleation_time < S.ymin] = 0
                f[nucleation_time > S.ymax] = 1

        return t, T, f


class Ferrite(PhaseTransformation):
    """
    Austenite to ferrite phase transformation
    """
    Q = 27500*4.184  # activation energy
    n1 = 0.41  # exponential factor 1
    n2 = 3  # exponential factor 2

    def initialize(self):
        self.comp_factor = self.alloy.FC  # composition factor for calculating transformation time
        self.Ts = self.alloy.Ae3  # transformation start temperature
        self.Tf = self.alloy.Bs  # transformation finish temperature


class Pearlite(PhaseTransformation):
    """
    Austenite to pearlite phase transformation
    """
    Q = 27500*4.184  # activation energy
    n1 = 0.32  # exponential factor 1
    n2 = 3  # exponential factor 2

    def initialize(self):
        self.comp_factor = self.alloy.PC  # composition factor for calculating transformation time
        self.Ts = self.alloy.Ae1  # transformation start temperature
        self.Tf = self.alloy.Bs  # transformation finish temperature


class Bainite(PhaseTransformation):
    """
    Austenite to bainite phase transformation
    """
    Q = 27500*4.184  # activation energy
    n1 = 0.29  # exponential factor 1
    n2 = 2  # exponential factor 2

    def initialize(self):
        self.comp_factor = self.alloy.BC  # composition factor for calculating transformation time
        self.Ts = self.alloy.Bs  # transformation start temperature
        self.Tf = self.alloy.Ms  # transformation finish temperature


class Martensite:
    """
    Athermal austenite to martensite transformation
    """

    def __init__(self, alloy):
        self.alloy = alloy
        self.Ts = self.alloy.Ms

        C, Mn, Si, Ni, Cr, Mo, Co = parse_comp(**self.alloy.w)
        self.alloy.alpha = 1e-3*(27.2 - (0.14*Mn + 0.21*Si + 0.11*Cr + 0.08*Ni + 0.05*Mo) - 19.8*(1-np.exp(-1.56*C)))

    def get_transformed_fraction(self, t, T, n=1000):
        """
        t, T : iterables
        Koistinen-Marburger equation
        """
        if len(t) > 3:
            # Fits T(t) by spline
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            # Uses linear interpolator
            t2T = interp1d(t, T)

        t = np.linspace(min(t), max(t), n)
        T = t2T(t)
        f = np.full(T.shape, 0, dtype=float)

        filtr = T < self.alloy.Ms
        if np.any(filtr):
            f[filtr] = 1 - np.exp(-self.alloy.alpha*(self.alloy.Ms - T[filtr]))
        return t, T, f


class TransformationDiagrams:
    colors_dict = dict(ferrite='#1f77b4', pearlite='#ff7f0e', bainite='#2ca02c',
                       martensite='#d62728', austenite='#9467bd')
    columns_label_dict = dict(t='Time (s)', T=u'Temperature (°C)',
                              ferrite='Ferrite', pearlite='Pearlite', bainite='Bainite',
                              martensite='Martensite', austenite='Austenite')

    def __init__(self, alloy):
        self.alloy = alloy

        self.ferrite = Ferrite(self.alloy)
        self.pearlite = Pearlite(self.alloy)
        self.bainite = Bainite(self.alloy)
        self.martensite = Martensite(self.alloy)

    def get_transformed_fraction(self, t, T, n=1000):
        """
        Get transformation curves for a given T(t) thermal cycle
        """
        t, T, f_ferr = self.ferrite.get_transformed_fraction(t, T, n)
        t, T, f_pear = self.pearlite.get_transformed_fraction(t, T, n)
        t, T, f_bain = self.bainite.get_transformed_fraction(t, T, n)
        t, T, f_mart = self.martensite.get_transformed_fraction(t, T, n)

        f_ferr_inc = np.diff(f_ferr, prepend=0)
        f_pear_inc = np.diff(f_pear, prepend=0)
        f_bain_inc = np.diff(f_bain, prepend=0)
        f_mart_inc = np.diff(f_mart, prepend=0)

        f_corr = pd.DataFrame(columns=['t', 'T', 'ferrite', 'pearlite', 'bainite', 'martensite', 'austenite'])
        f_corr['t'] = t
        f_corr['T'] = T
        f_corr.fillna(0, inplace=True)
        f_corr.loc[0, 'ferrite'] = f_ferr[0]
        f_corr.loc[0, 'pearlite'] = f_pear[0]
        f_corr.loc[0, 'bainite'] = f_bain[0]
        f_corr.loc[0, 'martensite'] = f_mart[0]
        f_corr.loc[0, 'austenite'] = 1. - f_ferr[0] - f_pear[0] - f_bain[0] - f_mart[0]

        def f1(i, x, y, z, w):
            if f_ferr[i] < 1:
                return f_corr.loc[i-1, 'ferrite'] + f_ferr_inc[i]*(1 - x - y - z - w)/(1 - f_ferr[i]) - x
            else:
                return f_corr.loc[i-1, 'ferrite'] + f_ferr_inc[i]*(1 - y - z - w) - x

        def f2(i, x, y, z, w):
            if f_pear[i] < 1:
                return f_corr.loc[i-1, 'pearlite'] + f_pear_inc[i]*(1 - x - y - z - w)/(1 - f_pear[i]) - y
            else:
                return f_corr.loc[i-1, 'pearlite'] + f_pear_inc[i]*(1 - x - z - w) - y

        def f3(i, x, y, z, w): return f_corr.loc[i-1, 'bainite'] + f_bain_inc[i]*(1 - x - y - w) - z

        def f4(i, x, y, z, w): return f_corr.loc[i-1, 'martensite'] + f_mart_inc[i]*(1 - x - y - z) - w

        for i in range(len(f_corr))[1:]:
            x0 = [f_corr.loc[i-1, 'ferrite'], f_corr.loc[i-1, 'pearlite'],
                  f_corr.loc[i-1, 'bainite'], f_corr.loc[i-1, 'martensite']]

            res = root(lambda x: [f1(i, *x), f2(i, *x), f3(i, *x), f4(i, *x)], x0=x0)

            f_corr.loc[i, 'ferrite'] = res.x[0]
            f_corr.loc[i, 'pearlite'] = res.x[1]
            f_corr.loc[i, 'bainite'] = res.x[2]
            f_corr.loc[i, 'martensite'] = res.x[3]
            f_corr.loc[i, 'austenite'] = 1. - res.x.sum()

        return f_corr.round(12)

    def draw_thermal_cycle(self, t, T, n=100, ax=None, **kwargs):
        """
        Draw thermal cycle (cooling curve) over AxesSubplot object
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        if len(t) > 3:
            # Fits T(t) by spline
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            # Uses linear interpolator
            t2T = interp1d(t, T)

        t = np.linspace(min(t), max(t), n)
        T = t2T(t)

        kw = dict(color='k', ls='--')
        kw.update(kwargs)

        return ax.plot(t, T, **kw)

    def TTT(self, fs=1e-2, ff=.99, ax=None, **kwargs):
        """
        Plot TTT diagram
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))
        else:
            fig = ax.get_figure()

        # Ferrite
        T = np.arange(self.alloy.Bs, self.alloy.Ae3)
        ts = self.ferrite.get_transformation_time(T, fs)  # start
        tf = self.ferrite.get_transformation_time(T, ff)  # finish
        ax.plot(ts, T, color=self.colors_dict['ferrite'], label='Ferrite {:g}%'.format(100*fs), **kwargs)
        ax.plot(tf, T, color=self.colors_dict['ferrite'], ls='--', label='Ferrite {:g}%'.format(100*ff), **kwargs)

        # Pearlite
        T = np.arange(self.alloy.Bs, self.alloy.Ae1)
        ts = self.pearlite.get_transformation_time(T, fs)
        tf = self.pearlite.get_transformation_time(T, ff)
        ax.plot(ts, T, color=self.colors_dict['pearlite'], label='Pearlite {:g}%'.format(100*fs), **kwargs)
        ax.plot(tf, T, color=self.colors_dict['pearlite'], ls='--', label='Pearlite {:g}%'.format(100*ff), **kwargs)

        # Bainite
        T = np.arange(self.alloy.Ms, self.alloy.Bs)
        ts = self.bainite.get_transformation_time(T, fs)
        tf = self.bainite.get_transformation_time(T, ff)
        ax.plot(ts, T, color=self.colors_dict['bainite'], label='Bainite {:g}%'.format(100*fs), **kwargs)
        ax.plot(tf, T, color=self.colors_dict['bainite'], ls='--', label='Bainite {:g}%'.format(100*ff), **kwargs)

        ax.axhline(self.alloy.Bs, color=self.colors_dict['bainite'], ls=':', label='Bs')
        ax.axhline(self.alloy.Ms, color=self.colors_dict['martensite'], label='Ms')

        ax.set_xscale('log')
        ax.set_ylim(25)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel(u'Temperature (°C)')
        ax.set_title(self.alloy.format_composition())
        ax.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, -.15))
        fig.tight_layout()

        return ax

    def CCT(self, Tini=900, fs=1e-2, ff=.99, cooling_rates=10**np.linspace(-4, 4, 320), ax=None, **kwargs):
        """
        Plot CCT diagram
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))
        else:
            fig = ax.get_figure()

        cooling_rates = np.array(cooling_rates)
        draw_cooling = kwargs.pop('draw_cooling', True)

        # Ferrite
        Ts = self.ferrite.get_transformation_temperature(Tini, self.alloy.Bs, cooling_rates, fs)  # start
        Tf = self.ferrite.get_transformation_temperature(Tini, self.alloy.Bs, cooling_rates, ff)  # finish
        ax.plot(Ts/cooling_rates, Ts, color=self.colors_dict['ferrite'], label='Ferrite {:g}%'.format(100*fs), **kwargs)
        ax.plot(Tf/cooling_rates, Tf, color=self.colors_dict['ferrite'],
                ls='--', label='Ferrite {:g}%'.format(100*ff), **kwargs)

        # Pearlite
        Ts = self.pearlite.get_transformation_temperature(Tini, self.alloy.Bs, cooling_rates, fs)
        Tf = self.pearlite.get_transformation_temperature(Tini, self.alloy.Bs, cooling_rates, ff)
        ax.plot(Ts/cooling_rates, Ts, color=self.colors_dict['pearlite'],
                label='Pearlite {:g}%'.format(100*fs), **kwargs)
        ax.plot(Tf/cooling_rates, Tf, color=self.colors_dict['pearlite'],
                ls='--', label='Pearlite {:g}%'.format(100*ff), **kwargs)

        # Bainite
        Ts = self.bainite.get_transformation_temperature(Tini, self.alloy.Ms, cooling_rates, fs)
        Tf = self.bainite.get_transformation_temperature(Tini, self.alloy.Ms, cooling_rates, ff)
        ax.plot(Ts/cooling_rates, Ts, color=self.colors_dict['bainite'], label='Bainite {:g}%'.format(100*fs), **kwargs)
        ax.plot(Tf/cooling_rates, Tf, color=self.colors_dict['bainite'],
                ls='--', label='Bainite {:g}%'.format(100*ff), **kwargs)

        ax.axhline(self.alloy.Bs, color=self.colors_dict['bainite'], ls=':', label='Bs')
        ax.axhline(self.alloy.Ms, color=self.colors_dict['martensite'], label='Ms')

        # Draw cooling curves
        if draw_cooling:
            for cooling_rate in cooling_rates[::10]:
                T = np.linspace(Tini, 25, 100)
                t = (Tini - T)/cooling_rate
                kw = dict(lw=.5)
                kw.update(kwargs)
                ax.plot(t, T, 'k:', **kw)

        ax.set_xscale('log')
        ax.set_ylim(25)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel(u'Temperature (°C)')
        ax.set_title(self.alloy.format_composition())
        ax.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, -.15))
        fig.tight_layout()

        return ax

    def plot_phase_fraction(self, t, T, n=1000, xaxis='t', ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        if len(t) > 3:
            # Fits T(t) by spline
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            # Uses linear interpolator
            t2T = interp1d(t, T)

        t = np.linspace(min(t), max(t), n)
        T = t2T(t)

        f = self.get_transformed_fraction(t, T, n)
        if f['ferrite'].max() > 0:
            ax.plot(f[xaxis], f['ferrite'], color=self.colors_dict['ferrite'], label='Ferrite')
        if f['pearlite'].max() > 0:
            ax.plot(f[xaxis], f['pearlite'], color=self.colors_dict['pearlite'], label='Pearlite')
        if f['bainite'].max() > 0:
            ax.plot(f[xaxis], f['bainite'], color=self.colors_dict['bainite'], label='Bainite')
        if f['martensite'].max() > 0:
            ax.plot(f[xaxis], f['martensite'], color=self.colors_dict['martensite'], label='Martensite')
        if f['austenite'].max() > 0:
            ax.plot(f[xaxis], f['austenite'], color=self.colors_dict['austenite'], label='Austenite')

        ax.set_xlabel(self.columns_label_dict[xaxis])
        ax.set_ylabel('Phase fraction')
        ax.legend()

        return ax


if __name__ == '__main__':
    # Defines alloy (grain size gs and composition)
    alloy = Alloy(gs=7, C=0.37, Mn=0.77, Si=0.15, Ni=0.04, Cr=0.98, Mo=0.21)

    # Initializes diagrams object
    diagrams = TransformationDiagrams(alloy)

    # Example 1: Plot TTT and CCT diagrams

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

    # Example 2: Plot CCT diagram and transformed fraction

    fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 6))
    fig2.subplots_adjust(wspace=.2)

    # Plot CCT diagram
    diagrams.CCT(Tini=1000, ax=ax3)

    # Chooses a thermal cycle (continuous cooling from 1000 to 0 oC in
    # 2000 s), draws cooling curve in the CCT and plots phase fraction
    t, T = [0, 2000], [1000, 0]
    # t, T = [0, 10000], [700, 700]
    diagrams.draw_thermal_cycle(t, T, ax=ax3)
    diagrams.plot_phase_fraction(t, T, xaxis='T', ax=ax4)

    fig2.suptitle(ax3.get_title())
    ax3.set_title('')
    ax4.set_title('')

    plt.show()
