import scipy.optimize as optim

class AlphaMinimizer:

    def __init__(self, maxiter=1000):
        self.minimizer = optim.minimize_scalar
        self.method = 'bounded'
        self.maxiter = maxiter

    def set_bounds(self, l_bound, u_bound):
        self.bounds = (l_bound, u_bound)
    
    def set_function(self, f):
        self.function = f

    def do_minimization(self, args):
        res = self.minimizer(
                self.function,
                bounds=self.bounds,
                args=args,
                method=self.method,
                options={'maxiter': self.maxiter}
                )
        if not res.success:
            return res.success, res.message
        else:
            return res.x, res.fun


