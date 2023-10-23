
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


import sys
sys.path.append('../')
import helpers
import parameters


pca = parameters.pca_default # nL/s
dR = parameters.dR # cm
n_sweep = 2


n_draw = 1001 # number of draws



with open('supporting/zeta_accuracy_tdiff.txt', 'w') as f:
    f.write('thold | t0 | tdiff | zeta | trend |\n')

count_t = 0
while count_t < n_draw:
    thold, t0, tdiff = helpers.draw_t()
    if min(thold, t0, tdiff) == tdiff:
            
        a = t0/thold
        b = t0/tdiff
        zeta = a*b

        _, peak_n2 = helpers.peak_ca(thold, t0, tdiff, n_sweep, pca, dR)

        if peak_n2 > 0:
            res = "saturating"
        else:
            res = "inverse"

        with open('supporting/zeta_accuracy_tdiff.txt', 'a') as f:
            f.write(f'{thold} | {t0} | {tdiff} | {zeta} | {res} |\n')
        f.close()

        count_t += 1


# array for store

with open('supporting/zeta_accuracy_thold.txt', 'w') as f:
    f.write('thold | t0 | tdiff | zeta | trend |\n')

count_t = 0
while count_t < n_draw:
    thold, t0, tdiff = helpers.draw_t()
    if min(thold, t0, tdiff) == thold:
            
        a = t0/thold
        b = t0/tdiff
        zeta = a*b

        _, peak_n2 = helpers.peak_ca(thold, t0, tdiff, n_sweep, pca, dR)

        if peak_n2 > 0:
            res = "saturating"
        else:
            res = "inverse"

        with open('supporting/zeta_accuracy_thold.txt', 'a') as f:
            f.write(f'{thold} | {t0} | {tdiff} | {zeta} | {res} |\n')
        f.close()

        count_t += 1






