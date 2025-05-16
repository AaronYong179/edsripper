# import zipfile

# def unzip_file(zip_filepath, extract_to_path):
#     with zipfile.ZipFile(zip_filepath, 'r') as zip_ref:
#         zip_ref.extractall(extract_to_path)

# # Example usage:
# zip_filepath = '2025-05-16_104732.eds'
# extract_to_path = '.'
# unzip_file(zip_filepath, extract_to_path)

import os

def convert_index_to_alphanum(well_index):
    """ assumes a 96 well plate format """
    column_number = well_index % 12 + 1
    row_alpha = chr(well_index // 12 + 65)
    return f"{row_alpha}{column_number}"

def read_data():
    fpath = os.path.join("apldbio", "sds", "analysis_result.txt")
    
    data = []
    raw_values = [
        [], # rn_values
        [], # delta_rn_values
    ]

    with open(fpath) as file:
        
        ## Ignore main header
        next(file) 

        ## Collect column headers
        # values are tab-separated. Trailing newline characters are stripped.  
        headers = list(map(lambda s : s.strip('\n'), next(file).split('\t')))
        for i, row in enumerate(file):
            if i % 3 == 0: # summary row
                data.append(list(map(lambda s : s.strip('\n'), row.split('\t'))))
            else:
                raw_values[i%3-1].append(list(map(float, row.split('\t')[1:])))

    print(raw_values[0])                    
    # return data

print(read_data())
# # print(convert_index_to_alphanum(93))
# d = list(map(float, read_data()[3][1:]))
# import matplotlib.pyplot as plt
# plt.plot(d)
# plt.show()