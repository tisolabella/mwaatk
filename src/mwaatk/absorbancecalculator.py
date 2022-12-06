import scipy.optimize as optim
import numpy as np
import matplotlib.pyplot as plt
import time

class ConvergenceError(Exception):
    pass


class AbsorbanceCalculator:

    def __init__(
            self, filter_object, optdict
            ):
        '''
        Initializes the calculator object. The object will have to perform the minimization of the function 
        F_1^p + F_2^p, where p = self.minimizer_power, which ideally should reach zero. 
        '''
        self.options = optdict
        self.filter = filter_object
        self.n_minimizations = self.options['n_min']
        self.tolerance = self.options['tolerance']
        self.minimizer_power = 2     #TODO- becomes an option in the config file
        algorithm = self.options['minimizer']

        if algorithm == 'dual_annealing':
            self.minimizer = optim.dual_annealing
            self.minimizer_options = {
                    'maxiter':self.options['max_iter'],
                    }
        elif algorithm == 'shgo':
            self.minimizer = optim.shgo
            self.minimizer_options = {
                    'options' : {'f_tol' : self.options['tolerance']},
                    'sampling_method' : 'sobol'
                    }
        elif algorithm == 'differential_evolution':
            self.minimizer = optim.differential_evolution
            self.minimizer_options = {
                    'maxiter' : self.options['max_iter'],
                    'popsize' : self.options['popsize'],
                    'disp' : False,
                    'updating' : self.options['updating'],
                    'workers' : self.options['workers'],
                }
        else: 
            print("RUNTIME ERROR: Double-check algorithm name")
        self.bounds_set = False

    def set_tolerance(self, tol):
        self.tolerance=tol

    def set_maxiter(self, mit):
        self.minimizer_options['maxiter'] = mit

    def set_bounds(self, omega_bounds, tau_bounds):
        self.omega_bounds = omega_bounds
        self.tau_bounds = tau_bounds
        lb = [omega_bounds[0], tau_bounds[0]]
        ub = [omega_bounds[1], tau_bounds[1]]
        self.bounds = optim.Bounds(lb, ub)
        self.bounds_set = True

    def set_filter(self, f):
        self.filter = f
    
    def do_minimization(self): 
        '''
        The method that carries out the minimization of the h_function.
        It returns a tuple containing the optimised parameter and the value
        of the minimized function (omega, tau, h_func).
        Since the differential evolution method is stochastic in nature,
        multiple minimizations are found to ascertain the stability of the
        solution. The number of these minimizations is controlled by the 
        instance variable self.n_minimizations, and the tolerance to accept
        the algorithm is self.tolerance
        '''
        if not self.bounds_set:
            self.set_bounds([0,1], self.options['tau_bounds'])
        results = []
        n_outliers = 0
        outliers_values = []
        omega = []
        tau = []
        for i in range(self.n_minimizations):
            res = self.minimizer(
                    self.h_function,
                    ((0,1), (0,20)),
                    **self.minimizer_options
                    )
            results.append(res) 
            if not res.success:
                if self.options['verbose']:
                    print("Minimization exited with status: {}, and reason : {}".format(res.success, res.message))
                return res.success
            elif res.fun >= 0 and res.fun < self.tolerance:
                omega.append(res.x[0])
                tau.append(res.x[1])
            else:
                n_outliers += 1 
                outliers_values.append(res.fun)
        if n_outliers > 0:
            if self.options['verbose'] and n_outliers < self.n_minimizations:
                print("Minimization succesful with warning: ignored {} minimization result(s) because the value of the function is further than tolerance from zero.\nThe outlier function values are {} ".format(n_outliers, outliers_values))
            elif self.options['verbose'] and n_outliers == self.n_minimizations:
                print("Minimization failed: ignored {} minimization result(s) because the value of the function is further than tolerance from zero.\nThe outlier function values are {} ".format(n_outliers, outliers_values))
        else:
            if self.options['verbose']:
                print("Minimization succesful. Ignored 0 minimization results. ")
        if omega != [] and tau != []:
            self.omega_min = np.mean(omega)
            self.tau_min = np.mean(tau)
            self.filter.converged = True
            return self.omega_min, self.tau_min 
        elif omega == [] or tau == []:
            return 'Error'

    def u_do_minimization(self): 
        if not self.bounds_set:
            self.set_bounds([0,1], self.options['tau_bounds'])
        results = []
        n_outliers = 0
        outliers_values = []
        omega = []
        tau = []
        if self.filter.converged:
            for i in range(self.n_minimizations):
                res = self.minimizer(
                        self.u_h_function,
                        ((0,1), (0,20)),
                        **self.minimizer_options
                        )
                results.append(res) 
                if not res.success:
                    self.filter.converged = False
                    if self.options['verbose']:
                        print("Minimization exited with status: {}, and reason : {}".format(res.success, res.message))
                    return res.success
                elif res.fun >= 0 and res.fun < self.tolerance:
                    omega.append(res.x[0])
                    tau.append(res.x[1])
                else:
                    n_outliers += 1 
                    outliers_values.append(res.fun)
            if n_outliers > 0:
                if self.options['verbose'] and n_outliers < self.n_minimizations:
                    print("Minimization succesful with warning: ignored {} minimization result(s) because the value of the function is further than tolerance from zero.\nThe outlier function values are {} ".format(n_outliers, outliers_values))
                elif self.options['verbose'] and n_outliers == self.n_minimizations:
                    print("Minimization failed: ignored {} minimization result(s) because the value of the function is further than tolerance from zero.\nThe outlier function values are {} ".format(n_outliers, outliers_values))
            else:
                if self.options['verbose']:
                    print("Minimization succesful. Ignored 0 minimization results. ")
            if omega != [] and tau != []:
                self.u_omega_min = np.mean(omega)
                self.u_tau_min = np.mean(tau)
                return self.u_omega_min, self.u_tau_min 
            elif omega == [] or tau == []:
                self.filter.converged = False
                return 'Error'           
        else:
            print(f"Uncertainty calculation for {self.filter.name} at {self.filter.wlength} nm not performed since the ABS minimization failed for this filter sample")
            return None, None

    def make_plots(self):
        om = np.linspace(0,1,5000)
        t  = np.linspace(0,30,5000)
        f1_om = self.f_1(om, self.tau_min)
        f1_t  = self.f_1(self.omega_min, t)
        f2_om = self.f_2(om, self.tau_min)
        f2_t  = self.f_2(self.omega_min, t)
        h_om = self.h_function((om, self.tau_min))
        h_t  = self.h_function((self.omega_min, t))
        #plots for fixed tau
        fig, (ax1, ax2)  = plt.subplots(2, 1, sharex=True)
        ax1.plot(om, f1_om, '-r', label='$F_1$')
        ax1.plot(om, f2_om, '-b', label='$F_2$')
        ax1.grid()
        ax1.legend()
        ax1.set_title(r'Fixed $\tau$ = {:.3f} for filter {} at wavelength {} nm'.format(self.tau_min, self.filter.name, self.filter.wlength))
        ax2.plot(om, h_om,'-k', label='$F_1^2 + F_2^2$')
        plt.xlabel('Single Scattering Albedo, $\omega$')
        plt.grid()
        plt.legend()
        plt.show()
        #plots for fixed omega 
        fig, (ax1, ax2)  = plt.subplots(2, 1, sharex=True)
        ax1.plot(t, f1_t, '-r', label='$F_1$')
        ax1.plot(t, f2_t, '-b', label='$F_2$')
        ax1.grid()
        ax1.legend()
        ax1.set_title(r'Fixed $\omega$ = {:.3f} for filter {} at wavelength {} nm'.format(self.omega_min, self.filter.name, self.filter.wlength))
        ax2.plot(t, h_t,'-k', label='$F_1^2 + F_2^2$')
        plt.xlabel('Optical thickness, $\\tau$')
        plt.grid()
        plt.legend()
        plt.show()

    def calc_abs(self):
        '''
        Calculate the sample absorbance with the parameters that minimize the h_function.
        '''
        try:
            if not self.filter.converged:
                raise ConvergenceError
            absorbance = (1 - self.omega_min) * self.tau_min
            absorbance = 100 * absorbance
            return absorbance
        except ConvergenceError:
            return None
   
    def u_calc_abs(self):
        '''
        Calculate the sample absorbance with the parameters that minimize the h_function.
        '''
        try:
            if not self.filter.converged:
                raise ConvergenceError
            absorbance = (1 - self.u_omega_min) * self.u_tau_min
            absorbance = 100 * absorbance
            return absorbance
        except ConvergenceError:
            return None


    def f_1(self, omega, tau, g=.75):
        '''
        f_1 is the first of the two functions which will be summed in quadrature and minimized.
        It is defined as the second member of Eq (1a) in Petzold 2003, minus the first member
        of the same equation.
        This  method f_1 when supplied the parameter g (set by default to the value 0.75).
        This method also needs to be passed an array of two variables. It will be only called within the
        minimization routine which will take care of the passing.
        From beta until F the parameters actually depend on the variables omega and tau,
        therefore it's necessary to calculate them each and every time the function is called
        '''
        beta = .5 * (1 - g - (4 / 25) * (1 - (np.abs(1 - 2 * g) / 8) - (7 / 8) * (1 - 2 * g) ** 2))
        beta_star = .5 * (1 - (g / 4) * (3 + g ** (3 + 2 * g ** 3)))
        a = 2 * (1 - omega * (1 - beta_star))
        b = 2 * omega * beta_star
        c = omega * beta
        d = omega * (1 - beta)
        p_1 = c - a * c - b * d
        p_2 = - a * d - b * c - d
        bigB = a ** 2 - b ** 2  #this is what Petzold 2003 calls B
        sqrtB = np.sqrt(bigB)
        T = np.exp(-tau)
        B_star = (b * (1 - T ** (2 * sqrtB))) / (sqrtB + a + (sqrtB - a) * T ** (2 * sqrtB))
        P_star = (1 / (2 * sqrtB) 
                * ((sqrtB - a + B_star * b) * T ** (-sqrtB) + (sqrtB + a - B_star * b) 
                * T ** sqrtB))
        B = (c - (p_1 / (1 + sqrtB)) - (c - (p_1 / (1 - sqrtB))) 
                * T ** (2 * sqrtB) - 2 * ((p_1 * sqrtB) / (1 - bigB)) 
                * T ** (1 + sqrtB)) / (sqrtB + a + (sqrtB - a) 
                * T ** (2 * sqrtB))
        F = (1 / (2 * sqrtB)) * ((d + B * b + (p_2 / (1 + sqrtB))) 
            * T ** (- sqrtB) - (d + B * b + (p_2 / (1 - sqrtB))) 
            * T ** sqrtB) + (p_2 / (1 - bigB)) * T

        # aggiungere commento di spiegazione
        B_m = self.filter.b_f_zero / (self.filter.p_f_zero + self.filter.b_f_zero)

        # The actual function
        f_1 = (T + F) / (1 - B_star * B_m) - self.filter.p_f / self.filter.p_f_zero
        return f_1


    def f_2(self, omega, tau, g=.75):
        '''
        f_2 is the second of the two functions which will be summed in quadrature and minimized.
        It is defined as the second member of Eq (1b) in Petzold 2003, minus the first member
        of the same equation.
        This method calculates f_2 when supplied the parameter g (set by default to the value 0.75).
        This method also needs to be passed an array of two variables. It will be only called within the
        minimization routine which will take care of the passing.
        From beta until F the parameters actually depend on the variables omega and tau,
        therefore it's necessary to calculate them each and every time the function is called
        '''
        beta = .5 * (1 - g - (4 / 25) * (1 - (np.abs(1 - 2 * g) / 8) - (7 / 8) * (1 - 2 * g) ** 2))
        beta_star = .5 * (1 - (g / 4) * (3 + g ** (3 + 2 * g ** 3)))
        a = 2 * (1 - omega * (1 - beta_star))
        b = 2 * omega * beta_star
        c = omega * beta
        d = omega * (1 - beta)
        p_1 = c - a * c - b * d
        p_2 = - a * d - b * c - d
        bigB = a ** 2 - b ** 2  #this is what Petzold 2003 calls B
        sqrtB = np.sqrt(bigB)
        T = np.exp(-tau)
        B_star = (b * (1 - T ** (2 * sqrtB))) / (sqrtB + a + (sqrtB - a) * T ** (2 * sqrtB))
        P_star = (1 / (2 * sqrtB) 
                * ((sqrtB - a + B_star * b) * T ** (-sqrtB) + (sqrtB + a - B_star * b) 
                * T ** sqrtB))
        B = (c - (p_1 / (1 + sqrtB)) - (c - (p_1 / (1 - sqrtB))) 
                * T ** (2 * sqrtB) - 2 * ((p_1 * sqrtB) / (1 - bigB)) 
                * T ** (1 + sqrtB)) / (sqrtB + a + (sqrtB - a) 
                * T ** (2 * sqrtB))
        F = (1 / (2 * sqrtB)) * ((d + B * b + (p_2 / (1 + sqrtB))) 
            * T ** (- sqrtB) - (d + B * b + (p_2 / (1 - sqrtB))) 
            * T ** sqrtB) + (p_2 / (1 - bigB)) * T

        # aggiungere commento di spiegazione
        B_m = self.filter.b_f_zero / (self.filter.p_f_zero + self.filter.b_f_zero)

        # The actual function
        f_2 = P_star * ((T + F) / (1 - B_star * B_m)) + B / B_m - self.filter.b_f / self.filter.b_f_zero
        return f_2
   

    def h_function(self, x, g=.75):
        '''
        This method conveniently wraps the function to minimize. It is simply the quadrature
        sum of f_1 and f_2.
        x is a 2D tuple containing omega and tau (necessary for making h callable by the minimizer).
        '''
        p = self.minimizer_power
        h_function = self.f_1(x[0], x[1], g) ** p + self.f_2(x[0], x[1], g) ** p
        return h_function

    def u_f_1(self, omega, tau, g=.75):
        '''
        f_1 is the first of the two functions which will be summed in quadrature and minimized.
        It is defined as the second member of Eq (1a) in Petzold 2003, minus the first member
        of the same equation.
        This  method f_1 when supplied the parameter g (set by default to the value 0.75).
        This method also needs to be passed an array of two variables. It will be only called within the
        minimization routine which will take care of the passing.
        From beta until F the parameters actually depend on the variables omega and tau,
        therefore it's necessary to calculate them each and every time the function is called
        '''
        beta = .5 * (1 - g - (4 / 25) * (1 - (np.abs(1 - 2 * g) / 8) - (7 / 8) * (1 - 2 * g) ** 2))
        beta_star = .5 * (1 - (g / 4) * (3 + g ** (3 + 2 * g ** 3)))
        a = 2 * (1 - omega * (1 - beta_star))
        b = 2 * omega * beta_star
        c = omega * beta
        d = omega * (1 - beta)
        p_1 = c - a * c - b * d
        p_2 = - a * d - b * c - d
        bigB = a ** 2 - b ** 2  #this is what Petzold 2003 calls B
        sqrtB = np.sqrt(bigB)
        T = np.exp(-tau)
        B_star = (b * (1 - T ** (2 * sqrtB))) / (sqrtB + a + (sqrtB - a) * T ** (2 * sqrtB))
        P_star = (1 / (2 * sqrtB) 
                * ((sqrtB - a + B_star * b) * T ** (-sqrtB) + (sqrtB + a - B_star * b) 
                * T ** sqrtB))
        B = (c - (p_1 / (1 + sqrtB)) - (c - (p_1 / (1 - sqrtB))) 
                * T ** (2 * sqrtB) - 2 * ((p_1 * sqrtB) / (1 - bigB)) 
                * T ** (1 + sqrtB)) / (sqrtB + a + (sqrtB - a) 
                * T ** (2 * sqrtB))
        F = (1 / (2 * sqrtB)) * ((d + B * b + (p_2 / (1 + sqrtB))) 
            * T ** (- sqrtB) - (d + B * b + (p_2 / (1 - sqrtB))) 
            * T ** sqrtB) + (p_2 / (1 - bigB)) * T

        # aggiungere commento di spiegazione
        B_m = self.filter.u_b_f_zero / (self.filter.u_p_f_zero + self.filter.u_b_f_zero)

        # The actual function
        f_1 = (T + F) / (1 - B_star * B_m) - self.filter.u_p_f / self.filter.u_p_f_zero
        return f_1


    def u_f_2(self, omega, tau, g=.75):
        '''
        f_2 is the second of the two functions which will be summed in quadrature and minimized.
        It is defined as the second member of Eq (1b) in Petzold 2003, minus the first member
        of the same equation.
        This method calculates f_2 when supplied the parameter g (set by default to the value 0.75).
        This method also needs to be passed an array of two variables. It will be only called within the
        minimization routine which will take care of the passing.
        From beta until F the parameters actually depend on the variables omega and tau,
        therefore it's necessary to calculate them each and every time the function is called
        '''
        beta = .5 * (1 - g - (4 / 25) * (1 - (np.abs(1 - 2 * g) / 8) - (7 / 8) * (1 - 2 * g) ** 2))
        beta_star = .5 * (1 - (g / 4) * (3 + g ** (3 + 2 * g ** 3)))
        a = 2 * (1 - omega * (1 - beta_star))
        b = 2 * omega * beta_star
        c = omega * beta
        d = omega * (1 - beta)
        p_1 = c - a * c - b * d
        p_2 = - a * d - b * c - d
        bigB = a ** 2 - b ** 2  #this is what Petzold 2003 calls B
        sqrtB = np.sqrt(bigB)
        T = np.exp(-tau)
        B_star = (b * (1 - T ** (2 * sqrtB))) / (sqrtB + a + (sqrtB - a) * T ** (2 * sqrtB))
        P_star = (1 / (2 * sqrtB) 
                * ((sqrtB - a + B_star * b) * T ** (-sqrtB) + (sqrtB + a - B_star * b) 
                * T ** sqrtB))
        B = (c - (p_1 / (1 + sqrtB)) - (c - (p_1 / (1 - sqrtB))) 
                * T ** (2 * sqrtB) - 2 * ((p_1 * sqrtB) / (1 - bigB)) 
                * T ** (1 + sqrtB)) / (sqrtB + a + (sqrtB - a) 
                * T ** (2 * sqrtB))
        F = (1 / (2 * sqrtB)) * ((d + B * b + (p_2 / (1 + sqrtB))) 
            * T ** (- sqrtB) - (d + B * b + (p_2 / (1 - sqrtB))) 
            * T ** sqrtB) + (p_2 / (1 - bigB)) * T

        # aggiungere commento di spiegazione
        B_m = self.filter.u_b_f_zero / (self.filter.u_p_f_zero + self.filter.u_b_f_zero)

        # The actual function
        f_2 = P_star * ((T + F) / (1 - B_star * B_m)) + B / B_m - self.filter.u_b_f / self.filter.u_b_f_zero
        return f_2
   

    def u_h_function(self, x, g=.75):
        '''
        This method conveniently wraps the function to minimize. It is simply the quadrature
        sum of f_1 and f_2.
        x is a 2D tuple containing omega and tau (necessary for making h callable by the minimizer).
        '''
        p = self.minimizer_power
        h_function = self.u_f_1(x[0], x[1], g) ** p + self.u_f_2(x[0], x[1], g) ** p
        return h_function


