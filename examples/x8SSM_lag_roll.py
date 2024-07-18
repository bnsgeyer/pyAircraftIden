import sys

sys.path.insert(0, '../')
from AircraftIden import FreqIdenSIMO, TransferFunctionFit
import math
import matplotlib.pyplot as plt
import pickle
import multiprocessing

import sympy as sp
from AircraftIden.StateSpaceIden import StateSpaceIdenSIMO, StateSpaceParamModel
import numpy as np
import csv
import sympy as sp

M = sp.Matrix([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, 1]])

g = 9.78

Yv, Yp, Ylat, wlg = sp.symbols('Yv Yp Ylat wlg')
Lv, Lp, Llat, wlag = sp.symbols('Lv Lp Llat wlag')


def callback(xk, state):
    print(xk)
    print(state)


def process_ssm(trims):
    with open('x8roll.csv', 'r') as f:
       reader = csv.reader(f)
       data = list(reader)
    arr = np.array(data)
    arr = np.array(data, dtype=float)

    time_seq_source = arr[:, 0]
    rout_source = arr[:, 1]
    gx_source = arr[:, 2]*math.pi / 180
    phi_source = arr[:, 3]*math.pi / 180
    ay_source = arr[:,4]
    vdot = arr[:,4]*3.282
    for i in range(len(ay_source)):
        vdot[i] = vdot[i] + phi_source[i]*32.2

    simo_iden = FreqIdenSIMO(time_seq_source, 1, 30, rout_source, vdot, gx_source, win_num=None)

    plt.rc("figure", figsize=(15,10))
    plt.figure("rout->vdot")
    simo_iden.plt_bode_plot(0)

    #plt.plot(time_seq_source, gx_source, color="red")
    plt.show()

    plt.figure("rout->p")
    simo_iden.plt_bode_plot(1)

    plt.show()

    ph0 = trims["ph0"]
    v0 = trims["v0"]
    F = sp.Matrix([[Yv, 0, 32.2, Ylat],
                   [Lv, 0, 0, Llat],
                   [0, 1, 0, 0],
                   [0, 0, 0, wlag]])
    G = sp.Matrix([[0], [0], [0], [wlg]])

    H0 = sp.Matrix([
        [0, 0, 0, 0],
        [0, 1, 0, 0]])
    H1 = sp.Matrix([
        [1, 0, 0, 0],
        [0, 0, 0, 0],
    ])
    syms = [Yv, Lv, Ylat, Llat, wlag, wlg]
    LatdynSSPM = StateSpaceParamModel(M, F, G, H0, H1, syms)
    con_str1 = ["A_3_3", "-B_3_0"]
    plt.rc('figure', figsize=(10.0, 5.0))

    freqres = simo_iden.get_freqres([0,1])
    for k in range(freqres.freq.__len__()):
        if (freqres.freq[k] > 10 or freqres.freq[k] < 2):
            freqres.coherens[0][k] = 0.0
        if (freqres.freq[k] > 30 or freqres.freq[k] < 1):
            freqres.coherens[1][k] = 0.0

    ssm_iden = StateSpaceIdenSIMO(freqres, accept_J=50,
                                  enable_debug_plot=False,
                                  y_names=["v", "p"], reg=0.1, iter_callback=callback, max_sample_times=1, con_str = con_str1)

# This allows the user to bound the random initial guess and the bound the minimize function when its called
    bnd = ([-0.84,  14, -0.8,  140, -40,30],[-0.44,15,-0.4,160,-30,40])
    initx0 = ([-0.64, 14.57, -0.58, 153, -35.3, 35.3])
    J, ssm = ssm_iden.estimate(LatdynSSPM, syms, constant_defines={}, bounds=bnd, initvals=initx0)

# This provides a random initial guess within +/- rand_init_max.  This does not bound the minimize function.
#    J, ssm = ssm_iden.estimate(LatdynSSPM, syms, constant_defines={}, rand_init_max=200)

    ssm.check_stable()
    ssm_iden.print_res()
    ssm_iden.draw_freq_res()

    plt.show()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    process_ssm({
        "ph0":0,
        "v0":0,
    })