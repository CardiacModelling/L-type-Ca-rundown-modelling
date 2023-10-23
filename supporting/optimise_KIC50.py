import copy
import numpy as np
import json
import matplotlib.pyplot as plt

import myokit
import pints  

# import project file
import helpers
import parameters

class EstimateIC50(pints.ForwardModel):
    def __init__(self, N_div):
        path = 'resources/zeng-ical-template.mmt'
        model = myokit.load_model(path)

        d_inf = model.get('L_type_Ca_current.d_infinity').pyfunc()(-90)
        f_inf = model.get('L_type_Ca_current.f_infinity').pyfunc()(-90)
        model.get('L_type_Ca_current.f').set_state_value(f_inf)
        model.get('L_type_Ca_current.d').set_state_value(d_inf)

        # ensure that this cell is at default value: 30 um
        model.set_value('calcium_dynamics.R', 30e-4)

        volt = model.get('membrane.V')
        volt.set_binding(None)
        B_arr = [0 for _ in range(N_div)]
        m = helpers.multi_compartment_model(model, N_div, B_arr)

        # create a copy to get initial conditions for buffer

        t0 = parameters.t0_default # milli seconds

        # array for store
        log_b = [f'calcium_dynamics.B{i}' for i in range(N_div)]

        m1 = copy.deepcopy(m)
        m1.get('calcium_dynamics.buffering').set_rhs(0)
        m1.get('calcium_dynamics.current').set_rhs(0)
        # strip ion-channel model
        m1.get('L_type_Ca_current.d').demote
        m1.get('L_type_Ca_current.d').set_rhs(0)
        m1.get('L_type_Ca_current.f').demote
        m1.get('L_type_Ca_current.f').set_rhs(1)

        for i in range(N_div):
            m1.get(f'calcium_dynamics.Ca{i}').demote
            m1.get(f'calcium_dynamics.Ca{i}').set_rhs(0)
            m1.get(f'calcium_dynamics.CaB{i}').demote
            m1.get(f'calcium_dynamics.CaB{i}').set_rhs(0)

        s = myokit.Simulation(m1)
        s.set_tolerance(abs_tol=parameters.tol_abs, rel_tol=parameters.tol_rel)

        d = s.run(t0 + 1, log = log_b, log_times = np.array([t0]))

        for b in log_b:
            m.get(b).set_state_value(d[b][0])

        # ----------------------------------


        self.sim = myokit.Simulation(m)
        self.sim.set_tolerance(abs_tol=parameters.tol_abs, rel_tol=parameters.tol_rel)


    def n_parameters(self):
        return 1

    def simulate(self, parameters, times):
        
        self.sim.reset()
        self.sim.set_constant('L_type_Ca_current.Km_Ca', parameters[0])

        try:
            
            self.sim.set_constant('membrane.V', -90)
            self.sim.run(10000, log = [])
            self.sim.set_constant('membrane.V', -0.001)
            d = self.sim.run(120, log = ['L_type_Ca_current.i_CaL'], log_times=times)
    
        except myokit.SimulationError:
            print('Error evaluating with parameters: ' + str(parameters))
            return np.nan * times

        return d['L_type_Ca_current.i_CaL']


initial_points = np.logspace(-7, -3, 3, endpoint=True)
boundary = pints.RectangularBoundaries(np.array([pow(10, -8)]), \
    np.array([pow(10, -2)]))

N_div0 = 600
initialmodel = EstimateIC50(N_div0)

times = np.arange(10000, 10120, 0.1)
syn_data = initialmodel.simulate(np.array([0.09e-5]), times)
noise = np.random.normal(0, 1, np.array(syn_data).shape)
syn_data +=noise


for N in np.arange(500, 700, 20):
    print(f'N is {N}')
    model_new = EstimateIC50(N)
    problem = pints.SingleOutputProblem(model_new, times, syn_data)
    error = pints.MeanSquaredError(problem)


    best_error = [None] * len(initial_points) 
    final_parameters = [None] * len(initial_points)
    time_taken = [None] * len(initial_points)

    for i in range(len(initial_points)):
        print(f"\n This is point {i + 1}")
        opt = pints.OptimisationController(error, initial_points[i], method=pints.NelderMead)
        opt.set_parallel(True)
        x2, f2 = opt.run()
        time_taken[i] = np.divide(opt.time(), 60)

        best_error[i] =f2
        final_parameters[i] = x2

        my_dict = {}
        my_dict['initial points'] = [x for x in initial_points]
        my_dict['final error score'] = best_error
        my_dict['final parameters'] = [x.tolist() for x in final_parameters[:i+1]]
        my_dict['time taken'] = time_taken

        with open(f'supporting/output/cdi_N{N}.json', 'w') as outfile:
                json.dump(my_dict, outfile, indent= 4)

