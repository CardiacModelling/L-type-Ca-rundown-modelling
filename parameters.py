
# Radius
R_min = 2e-4 #cm
R_max = 40.6e-4 #cm
R_default = 30e-4 # cm

# time before experiment
t0_min = 5000 # milli second
t0_max = 300000 # milli second
t0_default = 180000 # milli second

# thold
thold_min = 10000 # milli second
thold_max = 40000 # milli second
thold_default = 10000 # milli second

# PCa
pca_min = 0.005 # nL/s
pca_max = 0.0444 # 0.0457 # nL/s
pca_default = 0.04 # nL/s

# shell width
dR = 0.05e-4 # cm

# radius of the inner hemisphere
Rh = 1e-4/1.4142 # cm

# solver tolerance
tol_abs = 1e-9
tol_rel = 1e-9

# diffusion
DB = 2e-9 #cm2/ms
Bmax = 10 #mM

def cal_n_sweep(thold):
    "returns number of sweep and the multiplication factor"
    
    if thold == 10000:
        n_sweep = 20
        x = 10.12
    elif thold == 20000:
        n_sweep = 10
        x = 20.12
    elif thold == 40000:
        n_sweep = 5
        x = 40.12
    else:
        raise ValueError('thold not found')
    
    return n_sweep, x