## Original Ct mid derivation by Ian
## Modified by Aaron

import os, shutil, zipfile, csv
import logging, argparse
import numpy as np
from scipy.optimize import curve_fit

# setup logger
logging.basicConfig(
    level=logging.INFO, 
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

class Util:
    def _convert_index_to_alphanum(well_index:str) -> str:
        """ Converts a well index to its corresponding alphanumeric representation.

        For example, well index 8 = A9 and well index 20 = B9. This conversion method 
        assumes a 96 well format.
        """
        column_number = well_index % 12 + 1
        row_alpha = chr(well_index // 12 + 65)
        return f"{row_alpha}{column_number}"
    
    def _sigmoid(x, L, x0, k, b):
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
    def extract_eds_data(eds_path:str) -> str:
        """ Extracts only the raw analysis results from an .eds file. """

        logging.info("Preparing to extract analysis_result.txt file from .eds archive.")

        # check if master cached_results folder is present for saving the extracted results file
        if not os.path.isdir("cached_results"):
            os.mkdir("cached_results")

        # a quick check to ensure that the input path is actually an .eds file
        basename = os.path.basename(eds_path)
        filename, extension = os.path.splitext(basename)
        try:
            assert extension == ".eds"
        except:
            logging.error("Input file not of the expected .eds format.")
            raise Exception
        
        # if the file is already present (by name matching), then skip extraction
        OUTPUT_PATH = os.path.join("cached_results", f"{filename}_analysis_results.txt")
        if os.path.exists(OUTPUT_PATH):
            logging.info("Cached results already exists. Skipping extraction.")
            logging.warning("Make sure that all .eds file names are uniquely identifiable.")
            return OUTPUT_PATH

        # otherwise, extract only the analysis result file.
        with zipfile.ZipFile(eds_path, 'r') as ziphandle:
            RESULTS_PATH = "apldbio/sds/analysis_result.txt" # ensure that this matches the expected structure
            with ziphandle.open(RESULTS_PATH) as source, open(OUTPUT_PATH, "wb") as target:
                target.write(source.read())

        logging.info("Extraction complete")

        return OUTPUT_PATH # for further processing
    
    def default_parse_result_file(fhandle):

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
                well_key = Util._convert_index_to_alphanum(int(row[KEY_IDX]))
                data[well_key] = {}
                for j, header in enumerate(headers):
                    if j == KEY_IDX: continue
                    data[well_key][header] = row[j]
            else: # raw data rows
                assert well_key != None
                data[well_key][row[RAW_VALUE_KEY_IDX]] = list(map(float, row[1:]))

        return data
    
    def _fit_logito(rn_values):
        # The assumption here is that each Rn value corresponds to one cycle
        X = np.array(range(1, len(rn_values)+1))

        # A mandatory initial guess
        p0 = [
            max(rn_values) - min(rn_values), # guess L
            np.median(X), # guess x0
            1, # guess k
            min(rn_values) # guess b
        ]

        popt, _ = curve_fit(Util._sigmoid, X, rn_values, p0, method="lm")
        y_pred = Util._sigmoid(X, *popt) 
        return popt, y_pred
    
    def get_ct_values(data:dict) -> list:
        """ Obtain Ct values (both machine-returned and logistic fitted) and format as a table. 
        
        The output is a list of lists, where the first row is a header row. This makes it easy to write to a csv
        file.

        Args:
            data (dict) : parsed data from the input .eds file. 
        Returns:
            ct_values (list) : table formatted output     
        """

        # prepare header row 
        ct_values = [["Well", "Sample Name", "Detector", "Avg Ct", "Mid Ct"]]
        
        for key in data:
            popt, _ = Ripper._fit_logito(data[key]["Delta Rn values"])
            L, x0, k, b = popt # explicit unpacking for clarity
            ct_values.append(
                [key, data[key]["Sample Name"], data[key]["Detector"], float(data[key]["Avg Ct"]), float(x0)]
            )  
        return ct_values
    
    def write_rows(output_path:str, rows:list) -> None:
        """ Helper method that simply writes a list of lists to the specified output path"""
        with open(output_path, 'w', newline='') as wfile:
            csv.writer(wfile).writerows(rows)
        logging.info(f"Written to output file {output_path}")


def main():

    parser = argparse.ArgumentParser(
        prog="extract_eds", 
        description="Extracts raw qPCR from .eds files returned by Design and Analysis software.",
        add_help=True
    )

    parser.add_argument("-i", "--input", required=True, help="Path to the input .eds file")
    parser.add_argument("-o", "--output", required=True, help="Path to write the output file")
    parser.add_argument("-f", "--func", required=True, choices=["get_ct"], help="Function to run on input file")

    args = parser.parse_args()
    with open(Ripper.extract_eds_data(args.input)) as fhandle:
        data = Ripper.default_parse_result_file(fhandle)
        match args.func:
            case ("get_ct"):
                logging.info("Calling FUNCTION get_ct")
                Ripper.write_rows(args.output, Ripper.get_ct_values(data))

if __name__ == '__main__':
    main()