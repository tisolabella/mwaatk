class Filter:

    def __init__(self, p_f, p_f_zero, b_f, b_f_zero, name=None, wlength=None):
        ''' 
        As in the rest of the program, the names of the variables mirror the names they have 
        in the paper by Petzold (2003). For example, p_f is the fraction of radiation that passes
        the compound system aerosol layer + free filter matrix and so on.
        The area is in cm^2, while the volume is in liters.
        '''
        self.name = name
        self.p_f = p_f
        self.p_f_zero = p_f_zero
        self.b_f = b_f
        self.b_f_zero = b_f_zero
        self.wlength = wlength
        self.converged = False 

    def set_uncertainties(self, u_p_f, u_p_f_zero, u_b_f, u_b_f_zero):
        self.u_p_f = u_p_f
        self.u_p_f_zero = u_p_f_zero
        self.u_b_f = u_b_f
        self.u_b_f_zero = u_b_f_zero


    def print(self):
        print("Filter name:{}\nWavelength = {} nm \nPF / PF_0 = {:.4f} / {:.4f}\n BF / BF_0 = {:.4f} / {:.4f}".format(
            self.name, self.wlength, self.p_f, self.p_f_zero, self.b_f, self.b_f_zero))

