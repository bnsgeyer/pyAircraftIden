import sys

sys.path.insert(0, '../')

from AircraftIden import FreqIdenSIMO, TransferFunctionFit, TransferFunctionParamModel
from AircraftIden.TransferFunctionFit import plot_fitter
import matplotlib.pyplot as plt
import numpy as np
import csv
import sympy as sp
import math

def siso_freq_iden(win_num=32):
    with open('x8roll.csv', 'r') as f:
       reader = csv.reader(f)
       data = list(reader)
    arr = np.array(data)
    arr = np.array(data, dtype=float)
    #save_data_list = ["running_time", "yoke_pitch", "theta", "airspeed", "q", "aoa", "VVI", "alt"]
    #arr = np.load("../data/sweep_data_2017_11_16_11_47.npy")
    time_seq_source = arr[:, 0]
    rout_source = arr[:, 1]
    gx_source = arr[:, 2]*math.pi / 180
    simo_iden = FreqIdenSIMO(time_seq_source,0.5, 50, rout_source, gx_source, win_num=None)

    plt.rc("figure", figsize=(15,10))
    plt.figure("rout->p")
    simo_iden.plt_bode_plot(0)

    #plt.plot(time_seq_source, gx_source, color="red")

    plt.show()

    freq, H, gamma2, gxx, gxy, gyy = simo_iden.get_freq_iden(0)
    a, b, c, d, e, f, g, tau, s = sp.symbols("a b c d e f g tau s")
    num = f*s*s
    den = a*s*s*s*s + b*s*s*s + c*s*s + d*s + e
    tfpm = TransferFunctionParamModel(num, den, tau)
    fitter = TransferFunctionFit(freq, H, gamma2, tfpm, nw=20, iter_times=300, reg = 0.1)
    tf = fitter.estimate(2, 20, accept_J=50)
    plot_fitter(fitter, "$rout -> p$")

if __name__ == "__main__":
    siso_freq_iden(128)
