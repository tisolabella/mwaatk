from mwaatk.dataset import Dataset
from mwaatk.filtersample import Filter
import pandas as pd


class DatasetFromFile(Dataset):
    def __init__(self, fpath, optdict):
        self.options = optdict
        self.white_type = self.options['white_type']
        self.fpath = fpath

    def open_file(self):
        if self.fpath[-3:] == 'csv':
            self.dataset = self.open_csv()
        else:
            print("Invalid format")

    def open_csv(self):
        fulldata = pd.read_csv(self.fpath) 
        dataset = []
        for i, fname in fulldata['Name']:
            wl = fulldata['Wavelength'][i]
            p_f = fulldata['N1'][i]
            b_f = fulldata['N2'][i]
            p_f_zero = fulldata['B1'][i]
            b_f_zero = fulldata['B2'][i]
            dataset.append(Filter(p_f, p_f_zero, b_f, b_f_zero, wlength=wl,
                    name=fname, area=self.area, volume=self.volume)) 

        return dataset

    # All the other methods are defined in the mother class, Dataset

