import warnings
import json
from pprint import pprint
import csv
import numpy as np
from scipy.optimize import curve_fit
import scipy.odr as odr
import matplotlib.pyplot as plt

styledict = {0 : '.r', 1 : 'vb', 2 : '^g', 3 : 'sc', 4 : '3m', 5 : '<y', 6 : '>k', }


SCALE_GUESS = 10000.

def fit_function(x, scale, aae):
    return scale * x ** ( - aae)

def get_chisq(X, Y, unc_Y, scale, aae):
    chisq = 0
    for x, y, unc_y in zip(X, Y, unc_Y):
        chisq += ((y - fit_function(x, scale, aae)) ** 2) / (unc_y ** 2)
    return chisq

class AngstromCalculator:

    def __init__(self, data, opt):
        """
        The default constructor needs the data and the configuration. 
        Both can be either passed as a dictionary, if they come from 
        previous step in a larger program, or as .json files if this
        class is used as a standalone.
        
        Parameters:
            data: dict or string,
                either a dict containing the data in the MWAAT format,
                or an absolute path to a JSON file containing the data
                in MWAAT format;
            opt: dict or string,
                either a dict containing all the configuration options 
                or an absollute path to a JSON file containing all the
                configuration options.

        Returns:
            None
        """
        if type(opt) == dict:
            self.options = opt
        elif opt[-4:] == 'json':
            with open(opt, 'r') as f:
                self.options = json.load(f)
        if type(data) == dict:
            self.data = data
        elif data[-4:] == 'json':
            with open(data, 'r') as f:
                self.data = json.load(f)
        self.res = 5   #nm

    def writeout(self):
        if self.options['errors']:
            path = self.write_with_unc()
        else:
            path = self.write()
        return path

    def fit_angstrom(self):
        print("\n************************** \nCalculating Angstrom Absorption Exponents...")
        d = self.data   
        for fname in d.keys():
            ABS = []
            unc_ABS = []
            tmp_wl = []
            wlength = [float(x) for x in d[fname].keys() if 'ang' not in x]
            for w in wlength:
                if d[fname][str(int(w))]['ABS'] is not None:
                    ABS.append(float(d[fname][str(int(w))]['ABS']))
                    unc_ABS.append(float(d[fname][str(int(w))]['unc_ABS']))
                    tmp_wl.append(w)
            wlength = np.array(tmp_wl)
            unc_wl = self.options['wlength_res'] / np.sqrt(3) # Conversion to stat uncertainty
            unc_wl = np.array([unc_wl for w in wlength])
            ABS = np.array(ABS)
            unc_ABS = np.array(unc_ABS)
            # Set the fit bounds for the scale and aae
            fit_bounds = ([0, 0.5], [1e10, 5.]) 
            # Do fit using scipy.curve_fit
            fit_result = curve_fit(fit_function, wlength, ABS,
                    p0=(SCALE_GUESS, self.options['aae_guess']),
                    bounds=fit_bounds, sigma=unc_ABS)
            scale = fit_result[0][0]
            unc_scale = np.sqrt(fit_result[1][0][0])
            rel_unc_scale = unc_scale / scale * 100 
            aae = fit_result[0][1]
            unc_aae = np.sqrt(fit_result[1][1][1])
            rel_unc_aae = unc_aae / aae * 100 
            chisq = get_chisq(wlength, ABS, unc_ABS, scale, aae)
            red_chisq = chisq / 3 # There are 3 degrees of freedom
            if np.abs(aae) > 5. and np.abs(aae) < 10.:
                print(f"MWAAT-WARNING: the AAE  (fit result) for filter {fname} is {aae:.3f} +- {unc_aae:.3f}")
            elif np.abs(aae) > 10.:
                print(f"MWAAT-WARNING: the AAE (fit result) for filter {fname} is {aae:.3f} +/- {unc_aae:.3f}. Ignored. Try to change the guess in the options to obtain a better fit")
                d[fname]['angstrom'] = None
                d[fname]['u_angstrom'] = None
                continue
            d[fname]['angstrom'] = aae
            d[fname]['u_angstrom'] = unc_aae
            # Save the plot if required
            if self.options['abs_plot']:
                x = np.linspace(wlength[0], wlength[-1], 1000)
                y = fit_function(x, scale, aae)
                plt.subplots(2,1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
                ax1 = plt.subplot(2,1,1)
                ax1.errorbar(wlength, ABS, yerr=unc_ABS, xerr=unc_wl, fmt='none', 
                        ecolor='k', elinewidth=0.7, label='data')
                ax1.plot(x, y, '-b', linewidth=0.5, alpha=0.9, label='fit')
                ax1.annotate(r'Fit function: $C \lambda^{-AAE}$' + ' \n$AAE = $ {:.3f} $\pm$ {:.3f}\n'.format(aae, unc_aae) + r'$\chi^2_r = $' + f'{red_chisq:.3f}', xy=(0.65, 0.50), xycoords='axes fraction', bbox=dict(boxstyle="square,pad=0.5", fc='white', lw=2))
                ax1.set_title(f'Filter {fname}')
                ax1.set_ylabel(r'100 x ABS')
                ax1.grid()
                ax1.legend()
                # Calculate the residuals
                res_ABS = []
                for i, w in enumerate(wlength):
                    res_ABS.append( ABS[i] - fit_function(w, scale, aae) ) 
                ax2 = plt.subplot(2,1,2)
                ax2.errorbar(wlength, res_ABS, yerr=unc_ABS, xerr=unc_wl, fmt='none',
                        ecolor='k', elinewidth=0.7, label='residuals')
                ax2.set_xlabel('Wavelength [nm]')
                ax2.set_ylabel('Fit residuals')
                ax2.grid()
                pltpath = self.options['out_folder'] + 'MWAAT_plot_abs_' + str(fname) + '.png'
                plt.savefig(pltpath, dpi=500)
                plt.close()
    
    def write_with_unc(self):
        csvpath = self.options['out_folder'] + 'MWAAT_out_aae.csv'
        with open(csvpath, 'w') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(['Name', 'AAE'])
            for fn in self.data.keys():
                try:
                    ang = str(self.data[fn]['angstrom'])
                except KeyError:
                    ang = None
                row = [str(fn), ang]
                writer.writerow(row)
        return csvpath
    
    def write(self):
        csvpath = self.options['out_folder'] + 'MWAAT_out_aae.csv'
        with open(csvpath, 'r') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(['Name', 'AAE'])
            for fn in self.data.keys():
                ang = str(self.data[fn]['angstrom'])
                row = [str(fn), ang]
                writer.writerow(row)
        return csvpath


