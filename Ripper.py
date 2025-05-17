import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def sigmoid(x, L, x0, k, b): 
    return L / (1 + np.exp(-k*(x-x0))) + b

class Ripper:
    def __init__(self, fpath):
        self.fpath = fpath

    def _convert_index_to_alphanum(self, well_index):
        """ This conversion method assumes a 96 well format """
        column_number = well_index % 12 + 1
        row_alpha = chr(well_index // 12 + 65)
        return f"{row_alpha}{column_number}"
    

    def fit_logito(self, rn_values):
        X = np.array(range(1, len(rn_values)+1)) # assumed here that each Rn value corresponds to one cycle
        p0 = [max(rn_values) - min(rn_values), np.median(X), 1, min(rn_values)]
        popt, _ = curve_fit(sigmoid, X, rn_values, p0, method="lm")
        L, x0, k, b = popt
        y_pred = sigmoid(X, L, x0, k, b)
        plt.plot(rn_values)
        plt.plot(y_pred)
        plt.show()
        return x0
    
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
    
