import csv
import json
from copy import copy
from pprint import pprint
import numpy as np
import pandas as pd

from mwaatk.alphaminimizer import AlphaMinimizer

class Preprocessor:
    '''
    A class that contains all the routines to elaborate the raw data into passed radiation
    and backscattered radiation integrals.

    Attributes:
        options: dict
            a dictionary containing the configurations for the current analysis;
        pd_dis: float
            the distance between the scattering detecting photodiodes and the 
            sample;
        fwd_corr: float
            the corrective factor for the forward photodiode;
        bck_corr: float
            the integration factor for the backwards photodiodes;
        rho: float
            the filter surface roughness parameter;
        type: {'white', 'black'}
            the type of filter that the current Preprocessor instance is
            analysing;
        data: dict
            a multiple nested dictionary that contains all the raw data
            and the successive analysis results, up until the integrals;
        n_filters: int
            the number of filters in the sample;
        wave: sequence
            the wavelengths, preferably in nm, of the instrument;
        angles: sequence
            the angles in degrees of the photodiode positions;
    '''

    def __init__(self, typ, optdict):
        """
        Parameters:
            typ: {'white', 'black'}
                specifies the tipe of filters that this preprocessor instance will 
                elaborate;
            optdict: dict
                a dictionary containing the configurations for the current analysis.
        """
        self.options = optdict
        self.pd_dis = self.options['pd_distance']
        self.fwd_corr = self.options['fwd_corr']
        self.bck_corr = 2 * np.pi * self.pd_dis ** 2
        self.rho       = self.options['rho']   # The surface roughness parameter
        self.type = typ
        
    def raw_file_open(self, fpath, cols='A:D, H:J, N:P, T:V, Z:AB'):
        """
        Opens the file containing the raw data for analysis.

        Parameters:
            fpath: string,
                the file name containing the raw data formatted as in the 
                examples. Currently supports only Excel input, although
                a JSON input would be highly favoured;
            cols: string,
                the columns to read in the Excel input file. For info aboud
                formatting, see the documentation to pandas.read_excel.
        """
        if fpath[-3:] == 'xls' or fpath[-4:] == 'xlsx':
            self.data = self.open_excel(fpath, header=1, 
                    usecols=cols,
                    skiprows=0)
        elif fpath[-3:] == 'csv':
            self.data = self.open_csv(fpath) # altri args?
        elif fpath[-4:] == 'json':
            self.data = self.open_json(fpath) #altri args?

    def open_excel(self, fpath, header=None, usecols=None, skiprows=None):
        """
        Opens and read the data from the Excel input file.

        Parameters:
            fpath: string
                the file name containing the raw data formatted as in the 
                examples. 
            header: int
                the line in the xls(x) sheet to use as header;
            usecols: string,
                the columns to read;
            skiprows: int
                rows to skip at the beginning of the file.
        """
        rawdata = pd.read_excel(fpath, header=header, usecols=usecols, skiprows=skiprows)
        rawdata.head()
        # Sort names
        names = [x for x in rawdata['Name']]
        # Check whether the last filter is the first one 
        tmp = np.unique(np.array(names))
        if names[-1] == names[0] and len(tmp) > 1:
            names.append('_' + names[0])
        names =  np.unique(np.array(names))
        names = names.tolist()
        self.n_filters = len(names) 
        # Check whether wavelengths and angles are already set; if not, set the defaults
        # Create empty dictionary to fill tempdata
        try:
            self.wave
        except:
            self.set_lambdas(375, 407, 532, 635, 850)
        try:
            tmp_angle = [{} for x in self.angles]
        except:
            self.set_angles(0, 125, 165)
            tmp_angle = [{} for x in self.angles]
        # Create data dictionary
        tempdata = {f: 
                {w: 
                    {t: 
                        {
                            'raw': None, 
                            'mean': None, 
                            'gain': None, 
                            'dark': None,
                            'norm': None, 
                            'unc': None,
                            } 
                        for t in self.angles
                        } 
                    for w in self.wave
                    } 
                for f in names
                }
        # Fill the dictionary
        for k, values in rawdata.items():
            if k == 'Name':
                continue
            if k[-2:] == 'UV':
                wl = self.wave[0]
                theta = float(k[3:-4])
            if k[-2:] == ' B':
                wl = self.wave[1]
                theta = float(k[3:-3])
            if k[-2:] == ' G':
                wl = self.wave[2]
                theta = float(k[3:-3])
            if k[-2:] == ' R':
                wl = self.wave[3]
                theta = float(k[3:-3])
            if k[-2:] == 'IR':
                wl = self.wave[4]
                theta = float(k[3:-4])
            values = np.array(values)
            if self.options['white_type'] == 'full':
                for name in names:
                    indices = [j for j in rawdata.index if rawdata['Name'][j] == name] 
                    i_start = indices[0]
                    i_end = indices[-1] + 1
                    tempdata[name][wl][theta]['raw'] = copy(values[i_start:i_end])
                    #print(self.type, name, wl, theta, tempdata[name][wl][theta]['raw'])
                    #input()
            elif self.options['white_type'] == 'single' and self.type == 'black':
                for name in names:
                    indices = [j for j in rawdata.index if rawdata['Name'][j] == name] 
                    i_start = indices[0]
                    i_end = indices[-1] + 1
                    tempdata[name][wl][theta]['raw'] = copy(values[i_start:i_end])
                    #print(self.type, name, wl, theta, tempdata[name][wl][theta]['raw'])
                    #input()
            elif self.options['white_type'] == 'single' and self.type == 'white':
                # Do not split the data if there is only one white filter and
                # this is the white preprocessor.
                for name in names:
                    tempdata[name][wl][theta]['raw'] = copy(values)
                    #print(self.type, name, wl, theta, tempdata[name][wl][theta]['raw'])
                    #input()
        return tempdata 
        
    def open_csv(self, fpath):
        pass

    def open_json(self, fpath):
        pass

    def gd_file_open(self, fpath):
        """
        Opens the file containing the gains and dark current information for
        the photodiodes, formatted as in the examples.
        
        Parameters:
            fpath: string,
                the file path to the file containing the gains information.
                Currently only Excel is supported.
        """
        if fpath[-3:] == 'xls' or fpath[-4:] == 'xlsx':
            self.gd_open_excel(fpath)
        elif fpath[-3:] == 'csv':
            self.gd_open_csv(fpath) # altri args?
        elif fpath[-4:] == 'json':
            self.gd_open_json(fpath) #altri args?

    def gd_open_excel(self, fpath):
        gd = pd.read_excel(fpath, sheet_name=None)
        if self.type == 'white':
            ang = [int(x) for x in gd['wgd']['angle']]
            self.angles = ang
            lam = [int(x) for x in gd['wgd'].keys() if str(x) != 'angle']
            self.wave = lam
            empty = [{} for x in self.angles]
            gdict = {w: {th:x for (th,x) in zip(self.angles, copy(empty))} for w in self.wave}
            ddict = {w: {th:x for (th,x) in zip(self.angles, copy(empty))} for w in self.wave}
            for wl in self.wave:
                for th in self.angles:
                    #0 for 0°, 1 for 125°, 2 for 165°
                    gdict[wl][th] = gd['wgd'][wl][self.angles.index(th)] 
                    ddict[wl][th] = gd['wdd'][wl][self.angles.index(th)] 
        elif self.type == 'black':
            ang = [int(x) for x in gd['bgd']['angle']]
            self.angles = ang
            lam = [int(x) for x in gd['bgd'].keys() if str(x) != 'angle']
            self.wave = lam
            empty = [{} for x in self.angles]
            gdict = {w: {th:x for (th,x) in zip(self.angles, copy(empty))} for w in self.wave}
            ddict = {w: {th:x for (th,x) in zip(self.angles, copy(empty))} for w in self.wave}
            for wl in self.wave:
                for th in self.angles:
                    #0 for 0°, 1 for 125°, 2 for 165°
                    gdict[wl][th] = gd['bgd'][wl][self.angles.index(th)] 
                    ddict[wl][th] = gd['bdd'][wl][self.angles.index(th)] 
        self.set_gains(gdict)
        self.set_dark(ddict)

    def calculate_avg(self):
        """
        Updates self.data to contain another nested key , 'mean', which is
        the arithmetic average of the values in self.data[..][..][..]['raw'].
        """
        d = self.data
        for f_name in d.keys():
            for wlength in d[f_name].keys():
                for angle in d[f_name][wlength].keys():
                    dat = d[f_name][wlength][angle]['raw']
                    # Average the raw data
                    m = float(np.mean(dat))
                    d[f_name][wlength][angle]['mean'] = m 
        if self.options['errors']:
            self.calculate_uncertainty()
       
    def normalize_gain(self):
        """
        Normalizes the averaged data to account for different gains, according
        to the following method:
        let g_u be the highest gain among all lambdas and photodiodes; g_u will
        be the unit of gain to which every other value will be normalised.
        Then the normalized value V_ij for the j-th photodiode and the i-th 
        wavelength will be given by
        V_ij = (v_ij - d_ij) * (g_u / g_ij),
        where v_ij is the value prior to normalization, d_ij is the dark current
        for the wavelength-pdiode pair and g_ij related gain.
        The choice of g_u is arbitrary and could just as well be any of the gains, 
        or any number at all. A normalization to 1 might be formally more appropiate, but
        it is chosen to be the highest gain in order to streamline the method.

        The normalised value is then added to self.data with the key 'norm'.
        """
        g_u = []
        # Get a list of all gains:
        for wl in self.get_lambdas():
            for gain in self.get_gains(wl):
                g_u.append(gain)
        # Get the max gain:
        g_u = float(np.max(g_u))
        d = self.data
        for f_name in d.keys():
            for wlength in d[f_name].keys():
                for angle in d[f_name][wlength].keys():
                    # Dark subtraction and normalization:
                    v = d[f_name][wlength][angle]['mean']
                    g = d[f_name][wlength][angle]['gain']
                    dk = d[f_name][wlength][angle]['dark']
                    V = (v - dk) * (g_u / g)
                    d[f_name][wlength][angle]['norm'] = V

    def calculate_alpha(self):
        """
        This method calculates the fraction of diffuse radiation in the backwards
        hemisphere following the approach developed by Petzold 2003. 
        The alpha calculated for each filter, for each wavelength is added to the dict 
        at the level of the angle key. So the dict after this step will look like this:
        self.data = {
                'filter_name': {
                        lambda: {
                                angle: {
                                        'gain': g,
                                        ('raw': np.array([...]) #unless it was deleted
                                        'dark': d,
                                        'mean': m,
                                        'norm: n
                                        },
                                'alpha': a
                                }
                        }
                }                       
        """
        status = 1
        theta_1_rad = (self.angles[1] / 360.) * 2 * np.pi
        theta_2_rad = (self.angles[2] / 360.) * 2 * np.pi
        d = self.data
        m = AlphaMinimizer()
        m.set_bounds(0,1) # Alpha is a fraction
        m.set_function(self.f_alpha)
        for fname in d.keys():
            for wlength in d[fname].keys():
                x_1 = d[fname][wlength][125]['norm']
                x_2 = d[fname][wlength][165]['norm']
                args=(theta_1_rad, theta_2_rad, self.rho,
                        x_1, x_2) # The fixed arguments for each alpha minimization
                res = m.do_minimization(args)
                if not res[0]:
                    print('Minimizzazione non riuscita')
                    status = 0
                else:
                    alpha = float(res[0])
                    d[fname][wlength]['alpha'] = alpha
        if self.options['errors']:
            self.u_calculate_alpha()
        if status == 0:
            return False
        else:
            return True

    def u_calculate_alpha(self):
        """
        Calculates the diffuse fraction alpha using the boundary value for the raw data
        defined as the value plus its uncertainty.
        """
        status = 1
        theta_1_rad = (self.angles[1] / 360.) * 2 * np.pi
        theta_2_rad = (self.angles[2] / 360.) * 2 * np.pi
        d = self.data
        m = AlphaMinimizer()
        m.set_bounds(0,1) # Alpha is a fraction
        m.set_function(self.f_alpha)
        for fname in d.keys():
            for wlength in d[fname].keys():
                x_1 = d[fname][wlength][125]['norm'] + d[fname][wlength][125]['unc']
                x_2 = d[fname][wlength][165]['norm'] + d[fname][wlength][165]['unc']
                args=(theta_1_rad, theta_2_rad, self.rho,
                        x_1, x_2) # The fixed arguments for each alpha minimization
                res = m.do_minimization(args)
                if not res[0]:
                    print('Minimizzazione non riuscita')
                    status = 0
                else:
                    alpha = float(res[0])
                    d[fname][wlength]['u_alpha'] = alpha
        if status == 0:
            return False
        else:
            return True


    def calculate_piscattered(self):
        """ 
        Calculates the amount of radiation scattered at 180° knowing alpha and the radiation scattered at 165°.  
        """ 
        d = self.data
        theta_2_rad = (self.angles[2] / 360.) * 2 * np.pi
        for f in d.keys():
            for w in d[f].keys():
                alpha = d[f][w]['alpha']
                S_165 = d[f][w][165]['norm']
                cos = np.cos(theta_2_rad - np.pi)
                exponent = -.5 * (((theta_2_rad - np.pi) ** 2) / self.rho ** 2)
                exp = np.exp(exponent)
                S_180 = S_165 / (alpha * cos + (1 - alpha) * exp)
                d[f][w][180] = float(S_180) # Add to data dictionary

    def u_calculate_piscattered(self):
        """
        Same as calculate_piscattered but with the boundary value for raw data, defined
        as the value plus its uncertainty.
        """
        d = self.data
        theta_2_rad = (self.angles[2] / 360.) * 2 * np.pi
        for f in d.keys():
            for w in d[f].keys():
                alpha = d[f][w]['u_alpha']
                S_165 = d[f][w][165]['norm'] + d[f][w][165]['unc']
                cos = np.cos(theta_2_rad - np.pi)
                exponent = -.5 * (((theta_2_rad - np.pi) ** 2) / self.rho ** 2)
                exp = np.exp(exponent)
                S_180 = S_165 / (alpha * cos + (1 - alpha) * exp)
                d[f][w]['u_180'] = float(S_180) # Add to data dictionary


    def integrate(self):
        """
        Calculates the integrals for the passed radiation and the backscattered radiation, and saves
        them in the data dictionary.
        """
        d = self.data
        self.calculate_piscattered() # Radiation at 180°
        I_1 = 0.5  # Easy analytic integration
        theta = np.linspace((np.pi / 2), np.pi, 100000)
        exponent = - (((np.pi - theta) ** 2) / (2 * self.rho**2))
        int_2 = np.sin(theta) * np.exp(exponent)
        I_2 = np.trapz(int_2, x=theta)
        for f in d.keys():
            for w in d[f].keys():
                alpha = d[f][w]['alpha']
                S_0 = d[f][w][0]['norm']
                S_180 = d[f][w][180] / self.options['pd_area'] #current
                d[f][w]['p_f'] = S_0 / self.fwd_corr  # The passed radiation
                inte = self.bck_corr * (S_180 * (alpha * I_1 + (1 - alpha) * I_2)) # The backscattered radiation
                d[f][w]['b_f'] = inte
        if self.options['errors']:
            self.u_integrate()

    def u_integrate(self):
        """
        Calculates the integrals for the passed radiation and the backscattered radiation, and saves
        them in the data dictionary. In this case the calculations are performed on the value + uncertainty.
        """
        d = self.data
        self.u_calculate_piscattered() # Radiation at 180°
        I_1 = 0.5  # Easy analytic integration
        theta = np.linspace(np.pi / 2, np.pi, 10000)
        exponent = - (((np.pi - theta) ** 2) / (2 * self.rho**2))
        int_2 = np.sin(theta) * np.exp(exponent)
        I_2 = np.trapz(int_2, x=theta)
        for f in d.keys():
            for w in d[f].keys():
                alpha = d[f][w]['u_alpha']
                S_0 = d[f][w][0]['norm'] + d[f][w][0]['unc'] 
                S_180 = d[f][w]['u_180'] / self.options['pd_area'] #current
                d[f][w]['u_p_f'] = S_0 / self.fwd_corr  # The passed radiation
                d[f][w]['u_b_f'] = self.bck_corr * (S_180 * (alpha * I_1 + (1 - alpha) * I_2)) # The backscattered radiation
     
    def calculate_uncertainty(self):
        '''
        The uncertainty of each filter sample, at each wavelength
        for each photodiode is added to the dictionary under the key 'unc'.
        Then, the integrals will be recalculated repeating the whole process using 
        mean + unc for each data point.

        The methods starting with u_ are equivalent to the methods they refer to but they
        use the mean + unc instead of the mean for the calculations. Their definitions and
        implementations are immediately after each normal method.
        '''
        d = self.data
        for f_n in d.keys():
            for wl in d[f_n].keys():
                for ang in d[f_n][wl].keys():
                    if ang in self.angles:
                        # Upper dispersion uncertainty = max - mean
                        disp = float(np.max(d[f_n][wl][ang]['raw']) - np.mean(d[f_n][wl][ang]['raw']))
                        # Systematic uncertainty = dark current FOR NOW TAKEN OUT 
                        syst = float(d[f_n][wl][ang]['dark'])
                        uncertainty = disp #+ syst
                        d[f_n][wl][ang]['unc'] = uncertainty

    def csv_write(self, fpath):
        with open(fpath, 'w') as fwrite:
            writer = csv.writer(fwrite)
            header = ['Name', 'Type', 'Wavelength [nm]', 'Alpha', 'P_f', 'B_f']
            writer.writerow(header)
            for fname in self.data.keys():
                for wlength in self.data[fname].keys():
                    alpha, p_f, u_p_f, b_f, u_b_f = self.get_result(fname, wlength)
                    writer.writerow([fname, self.type, wlength, alpha, p_f, b_f])

    def json_write(self, fpath):
        names = self.get_names()
        wavelengths = self.get_wl()
        data_to_write = {n : {wl : {
                                    'alpha': None,
                                    'p_f': None,
                                    'u_p_f' : None,
                                    'b_f': None,
                                    'u_b_f' : None,
                                    'type': None
                                    }
                                for wl in wavelengths
                                }
                            for n in names
                            }
        for n in self.data.keys():
            for w in self.data[n].keys():
                    alpha, p_f, u_p_f, b_f, u_b_f = self.get_result(n, w)
                    data_to_write[n][w]['alpha'] = alpha
                    data_to_write[n][w]['p_f'] = p_f
                    data_to_write[n][w]['u_p_f'] = u_p_f
                    data_to_write[n][w]['b_f'] = b_f
                    data_to_write[n][w]['u_b_f'] = u_b_f
                    data_to_write[n][w]['type'] = self.type
        with open(fpath, 'w') as fwrite:
            json.dump(data_to_write, fwrite, indent=4)


    def s_alpha(self, alpha, theta=0.5, rho=0.5):
        """
        This method is used to calculate S(alpha | theta) as detailed in Petzold 2003.
        S is the functional form of the backscattered radiation assuming both diffuse 
        and Gaussian behaviour, in proportions dictated by the parameter alpha. 
        In this implementation, since alpha is the variable to to change in order to 
        find the zero later on, theta is the parameter and alpha the passed variable.
        
        """
        cos_t = np.cos(theta - np.pi)
        exp_t = np.exp(-0.5 * ((theta - np.pi) ** 2) / (rho ** 2))
        s_a = alpha * cos_t + (1 - alpha) * exp_t
        return s_a

    def r_alpha(self, alpha, theta_1=2.18, theta_2=2.88, rho=0.5):
        """
        The S-function ratio between S(alpha | theta = 125°) and S(alpha| theta = 165°)

        kwargs have to be theta_1=float, theta_2=float, rho=float
        """
        s_a_1 = self.s_alpha(alpha, theta=theta_1, rho=rho)
        s_a_2 = self.s_alpha(alpha, theta=theta_2, rho=rho)
        r_a   = s_a_1 / s_a_2
        return r_a

    def f_alpha(self, alpha, theta_1=2.18, theta_2=2.88, rho=0.5, val_1=5., val_2=8. ):
        """
        The difference between r_alpha and the measured ratio of the normalised photodiode
        readings at angles 125° and 165°.
        The method returns this difference squared to allow root finding through minimisation.

        kwargs have to be theta_1=float, theta_2=float, rho=float, val_1=float, val_2=float
        where val_1 and val_2 are the normalised photodiode readings at 125° and 165°,
        respectively.
        """
        r_a = self.r_alpha(alpha, theta_1=theta_1, theta_2=theta_2, rho=rho)
        r_x = val_1 / val_2
        f_a = r_x - r_a
        return f_a ** 2

    ### Setters #########
    def set_lambdas(self, *args):
        """
        this function needs to be passed all the wavelengths at which the
        instrument operates. For example self.set_lambda(420, 546, 635)
        """
        self.wave = args
    
    def set_angles(self, *args):
        """
        this function needs to be passed all the wavelengths at which the
        instrument operates. For example self.set_angles(420, 546, 635)
        """
        self.angles = args

    def set_gains(self, g_dict):
        """
        g_dict = {wl: 
                    {angle1: gain, ... 
                    }
                }
        """
        for wl in self.wave:
            for th in self.angles:
                for n in self.data.keys():
                    g = g_dict[wl][th]
                    self.data[n][wl][th]['gain'] = g

    def set_dark(self, d_dict):
        """
        d_dict = {wl: 
                    {angle1: dark, ... 
                    }
                }
        """
        for wl in self.wave:
            for th in self.angles:
                for n in self.data.keys():
                    d = d_dict[wl][th]
                    self.data[n][wl][th]['dark'] = d

    def set_theta_1(self, theta_1):
        """If one needed to change the first backward angle
        """
        self.angles[1] = theta_1

    def set_theta_2(self, theta_2):
        """If one needed to change the second backward angle
        """
        self.angles[2] = theta_2
     
    ### Getters #######
    def get_names(self):
        n = [x for x in self.data.keys()]
        return n

    def get_wl(self):
        w = []
        for n in self.data.keys():
            for wl in self.data[n].keys():
                if wl not in w:
                    w.append(wl)
        return w

    def get_lambdas(self):
        """
        Returns a tuple containing the wavelengths at which the instrument operates.
        """
        lam = []
        try:
            name = list(self.data.keys())[0]
            for w in self.data[name].keys():
                lam.append(w)
            return lam
        except:
            print("ERROR: data not yet read from file, or problem reading data")

    def get_gains(self, wlength):
        """
        Returns a list containing the photodiode gains for the specified lambda.
        (assuming equal gains for all filters at a given wavelength).
        """
        gain = []
        try:
            name = list(self.data.keys())[0]
            for a in self.data[name][wlength].keys():
                gain.append(self.data[name][wlength][a]['gain'])
            return gain
        except:
            print("ERROR: data not yet read from file, or problem reading data")

    def get_mean(self, fname, wlength):
        """
        Returns a list containing the photodiode mean for the specified lambda 
        and filter name
        """
        mean = []
        try:
            name = fname
            for a in self.data[name][wlength].keys():
                mean.append(self.data[name][wlength][a]['mean'])
            return mean
        except:
            print("ERROR: data not yet read from file, problem reading data, or \
                    data not yet averaged using calculate_avg")

    def get_all_data(self):
        """
        Returns the entire data dictionary.
        """
        return self.data

    def get_one_data(self, f, wl):
        """
        Returns the results from all the preprocessing stages for a 
        specified filter sample at a given wavelength.
        """
        return self.data[f][wl]

    def get_result(self, f, wl):
        """
        Returns the results of the analysis for the specified filter sample
        and wavelength.
        """
        d = self.data
        return d[f][wl]['alpha'], d[f][wl]['p_f'], d[f][wl]['u_p_f'], d[f][wl]['b_f'], d[f][wl]['u_b_f']



