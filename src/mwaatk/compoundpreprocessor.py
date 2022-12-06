from mwaatk.preprocessor import Preprocessor

import json
import numpy as np

class CompoundPreprocessor:
    """
    Abstraction level for the single preprocessors. It is necessary in order to normalize 
    the gains simultaneously, which is fundamental for consistency.
    """
    def __init__(self, optdict):
        self.options = optdict
        self.w_preprocessor = Preprocessor('white', optdict)
        self.b_preprocessor = Preprocessor('black', optdict)

    def set_lambdas(self, *args):
        self.w_preprocessor.set_lambdas(*args)
        self.b_preprocessor.set_lambdas(*args)

    def set_angles(self, *args):
        self.w_preprocessor.set_angles(*args)
        self.b_preprocessor.set_angles(*args)

    def raw_file_open(self):
        w_fpath = self.options['wfile']
        b_fpath = self.options['bfile']
        self.w_preprocessor.raw_file_open(w_fpath)
        self.b_preprocessor.raw_file_open(b_fpath)

    def calculate_avg(self):
        print("\n************************** \nStarting analysis...\n")
        self.w_preprocessor.calculate_avg()
        self.b_preprocessor.calculate_avg()

    def set_gains_dark(self):
        gains_dark_file = self.options['gfile']
        self.w_preprocessor.gd_file_open(gains_dark_file)
        self.b_preprocessor.gd_file_open(gains_dark_file)
    
    def calculate_uncertainty(self):
        self.w_preprocessor.calculate_uncertainty()
        self.b_preprocessor.calculate_uncertainty()

    def confirm_start(self):
        print(f"\nOpening \"{self.options['bfile']}\" for black data and \"{self.options['wfile']}\" for white data. White type is \"{self.options['white_type']}\".")
        print(f"Opening \"{self.options['gfile']}\" for gains and dark data")
        print(f"Tau bounds are set to {self.options['tau_bounds']}")
        if self.options['errors']:
            print(f"The uncertainty calculation will be performed.")
        print(f'The output files and plots will be saved in the folder {self.options["out_folder"]}.')
        print("Press ENTER to start the analysis, or CTRL-C to abort")

    def normalize_gains(self):
        """
        Takes care of normalizing the gains across both black and white.
        """
        g_u = []
        # Get a list of white gains:
        for wl in self.w_preprocessor.get_lambdas():
            for gain in self.w_preprocessor.get_gains(wl):
                g_u.append(gain)
        # Append all the black gains:
        for wl in self.b_preprocessor.get_lambdas():
            for gain in self.b_preprocessor.get_gains(wl):
                g_u.append(gain)
        # Get the max gain:
        g_u = np.max(g_u)
        # Normalize the white gains
        d = self.w_preprocessor.data
        for f_name in d.keys():
            for wlength in d[f_name].keys():
                for angle in d[f_name][wlength].keys():
                    # Dark subtraction and normalization:
                    v = d[f_name][wlength][angle]['mean']
                    g = d[f_name][wlength][angle]['gain']
                    dk = d[f_name][wlength][angle]['dark']
                    V = (v - dk) * (g_u / g)
                    d[f_name][wlength][angle]['norm'] = V
                    if self.options['errors']:
                        # Normalize the uncertainties too
                        u_v = d[f_name][wlength][angle]['unc']
                        g = d[f_name][wlength][angle]['gain']
                        u_V = (u_v) * (g_u / g)
                        d[f_name][wlength][angle]['unc'] = u_V
        # Normalize the black gains
        d = self.b_preprocessor.data
        for f_name in d.keys():
            for wlength in d[f_name].keys():
                for angle in d[f_name][wlength].keys():
                    # Dark subtraction and normalization:
                    v = d[f_name][wlength][angle]['mean']
                    g = d[f_name][wlength][angle]['gain']
                    dk = d[f_name][wlength][angle]['dark']
                    V = (v - dk) * (g_u / g)
                    d[f_name][wlength][angle]['norm'] = V
                    if self.options['errors']:
                        # Normalize the uncertainties too
                        u_v = d[f_name][wlength][angle]['unc']
                        g = d[f_name][wlength][angle]['gain']
                        u_V = (u_v) * (g_u / g)
                        d[f_name][wlength][angle]['unc'] = u_V


    def calculate_alpha(self):
        self.w_preprocessor.calculate_alpha()
        self.b_preprocessor.calculate_alpha()

    def integrate(self):
        self.w_preprocessor.integrate()
        self.b_preprocessor.integrate()
        print("\n**************************\nDone preprocessing.")

    def get_all_data(self):
        return self.b_preprocessor.get_all_data(), self.w_preprocessor.get_all_data()


    def write(self):
        bpath_json = self.options['out_folder'] + 'MWAAT_out_pre_b.json'
        wpath_json = self.options['out_folder'] + 'MWAAT_out_pre_w.json'
        bpath_csv = self.options['out_folder'] + 'MWAAT_out_pre_b.csv'
        wpath_csv = self.options['out_folder'] + 'MWAAT_out_pre_w.csv'
        self.b_preprocessor.csv_write(bpath_csv)
        self.w_preprocessor.csv_write(wpath_csv)
        self.b_preprocessor.json_write(bpath_json)
        self.w_preprocessor.json_write(wpath_json)
        return bpath_json, wpath_json

