## Original author: Ian Cheong
## Modified by: Aaron Yong

import numpy as np
import os, shutil, zipfile
from scipy.optimize import curve_fit

def sigmoid(x, L, x0, k, b):
    """
    The logistic function is defined as follows:
    f(x) = [ L / 1 + e^{-k(x - x_0)} ] + b

    where
        x is the cycle number (i.e, for a standard run, 1 <= x <= 40)
        L is the maximum signal at plateau
        x_0 is the cycle number at the midpoint of the curve
        k is the growth rate
        b is the minimum signal (baseline)

    Generally, the key parameters of interest should be the following:
        b, the minimum signal
        x_0, the cycle number at the midpoint of the curve, 
        y-gain, defined as L - b (i.e., max - min)
        PCR efficiency, defined as e^k
    """
    return L / (1 + np.exp(-k*(x-x0))) + b

class Ripper:
    def __init__(self, eds_path):

        ## Handle extraction first
        if not os.path.isdir("_temp"):
            os.mkdir("_temp")
        shutil.rmtree("_temp")

        if not os.path.isdir("results"):
            os.mkdir("results")
        
        basename = os.path.basename(eds_path)
        filename, extension = os.path.splitext(basename)
        assert extension == ".eds"

        with zipfile.ZipFile(eds_path, 'r') as ziphandle:
            ziphandle.extractall("_temp")
        
        analysis_results_path = os.path.join(
            ".", "_temp", "apldbio", "sds", "analysis_result.txt"
        )
        self.fpath = os.path.join(".", "results", f"{filename}_analysis_result.txt")
        shutil.copy(analysis_results_path, self.fpath)

        ## cleanup
        shutil.rmtree("_temp")


    def _convert_index_to_alphanum(self, well_index):
        """ This conversion method assumes a 96 well format """
        column_number = well_index % 12 + 1
        row_alpha = chr(well_index // 12 + 65)
        return f"{row_alpha}{column_number}"
    
    def fit_logito(self, rn_values):
        """ Fits a logistic (more specifically, a sigmoidal) curve to the raw qPCR data. 
        
        """

        # The assumption here is that each Rn value corresponds to one cycle
        X = np.array(range(1, len(rn_values)+1))

        # A mandatory initial guess
        p0 = [
            max(rn_values) - min(rn_values), # guess L
            np.median(X), # guess x0
            1, # guess k
            min(rn_values) # guess b
        ]

        popt, _ = curve_fit(sigmoid, X, rn_values, p0, method="lm")
        y_pred = sigmoid(X, *popt) 
        return popt, y_pred
    
    def _default_parse_result_file(self, fhandle):

        ## CONSTANTS
        KEY_IDX = 0 # well index number
        ITEM_ROWS = 3 # three rows per well (summary, Rn, and delta Rn)
        RAW_VALUE_KEY_IDX = 0
        
        ## Ignore main header
        next(fhandle)

        ## Collect column headers
        headers = list(map(lambda s : s.strip('\n'), next(fhandle).split('\t')))
        
        ## Collect data
        data = {}
        well_key = None
        for i, row in enumerate(fhandle):
            row = list(map(lambda s : s.strip('\n'), row.split('\t')))
            if i % ITEM_ROWS == 0: # summary row
                well_key = self._convert_index_to_alphanum(int(row[KEY_IDX]))
                data[well_key] = {}
                for j, header in enumerate(headers):
                    if j == KEY_IDX: continue
                    data[well_key][header] = row[j]
            else: # raw data rows
                assert well_key != None
                data[well_key][row[RAW_VALUE_KEY_IDX]] = list(map(float, row[1:]))

        return data
    
    def get_ct_values(self):
        
        with open(self.fpath) as fhandle:
            data = self._default_parse_result_file(fhandle)
        
        ct_values = {}
        for key in data:
            
            popt, _ = self.fit_logito(data[key]["Delta Rn values"])
            L, x0, k, b = popt # explicit unpacking for clarity
            
            ct_values[key] = { 
                "Avg Ct" : float(data[key]["Avg Ct"]),
                "Mid Ct" : float(x0)
            }

        print(ct_values)
        

    
