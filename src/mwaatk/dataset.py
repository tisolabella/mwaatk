import json
from mwaatk.filtersample import Filter

class Dataset:
    def __init__(self, black_dict, white_dict, optdict):
        '''
        The preprocessed data is passed to the Dataset object as dictionaries
        (obtained by two preprocessor objects, one for white and one for black)
        '''
        self.options = optdict
        self.white_type = self.options['white_type']

        if type(white_dict) == dict:
            self.white_dict = white_dict
        elif white_dict[-4:] == 'json':
            with open(white_dict, 'r') as fw:
                self.white_dict = json.load(fw)
        if type(black_dict) == dict:
            self.black_dict = white_dict
        elif black_dict[-4:] == 'json':
            with open(black_dict, 'r') as fb:
                self.black_dict = json.load(fb)

    def create_set(self):
        if self.white_type == 'single':
            self.dataset = self.single_white()
        elif self.white_type == 'full':
            self.dataset = self.full_white()

    def single_white(self):
        b_data = self.black_dict
        w_data = self.white_dict
        dataset = []
        # Get the name of the white filter
        w_fname = str(list(w_data.keys())[0]) 
        # Create all the filter objects
        for b_fname, b_val in b_data.items():
            for b_wl in b_val.keys():
                p_f = b_val[b_wl]['p_f']
                p_f_zero = w_data[w_fname][b_wl]['p_f']
                b_f = b_val[b_wl]['b_f']
                b_f_zero = w_data[w_fname][b_wl]['b_f']
                dataset.append(Filter(p_f, p_f_zero, b_f, b_f_zero, wlength=b_wl,
                        name=b_fname)) 
                if self.options['errors']:
                    u_p_f = b_val[b_wl]['u_p_f']
                    u_p_f_zero = w_data[w_fname][b_wl]['u_p_f']
                    u_b_f = b_val[b_wl]['u_b_f']
                    u_b_f_zero = w_data[w_fname][b_wl]['u_b_f']
                    dataset[-1].set_uncertainties(u_p_f, u_p_f_zero, u_b_f, u_b_f_zero)
        return dataset 

    def full_white(self):
        b_data = self.black_dict
        w_data = self.white_dict
        # Initialize empty sequence to hold the Filter objects
        dataset = []
        # Create the Filter objects
        for b_fname, b_val in b_data.items():
            for b_wl in b_val.keys():
                p_f = b_val[b_wl]['p_f']
                p_f_zero = w_data[b_fname][b_wl]['p_f']
                b_f = b_val[b_wl]['b_f']
                b_f_zero = w_data[b_fname][b_wl]['b_f']
                u_p_f = b_val[b_wl]['u_p_f']
                u_p_f_zero = w_data[b_fname][b_wl]['u_p_f']
                u_b_f = b_val[b_wl]['u_b_f']
                u_b_f_zero = w_data[b_fname][b_wl]['u_b_f']
                dataset.append(Filter(p_f, p_f_zero, b_f, b_f_zero, wlength=b_wl,
                        name=b_fname)) 
                dataset[-1].set_uncertainties(u_p_f, u_p_f_zero, u_b_f, u_b_f_zero)
        return dataset #a sequence of Filter object

    ### Setters #####

    ### Getters ####

    def get_full_set(self):
        """ 
        Return the entire dataset.
        """
        return self.dataset

    def get_filter(self, fname):
        """
        Return the data recorded in the filter with the given name.
        """
        flist = []
        for f in self.dataset:
            if f.name == fname:
                flist.append(f)
        return flist


