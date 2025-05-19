import os
import matplotlib.pyplot as plt
import numpy as np
from Ripper import Ripper

fpath = os.path.join("apldbio", "sds", "analysis_result.txt")
ripper = Ripper("2025-05-16_104732.eds")

# ripper.get_ct_values()

# with open(fpath) as fhandle:
#     data = ripper._default_parse_result_file(fhandle)
#     key = "B9"
#     x0 = ripper.fit_logito(data[key]["Delta Rn values"])
#     # plt.plot(data[key]["Delta Rn values"])
#     # plt.axvline(x0)
#     # plt.show()

# def sigmoid(x, L, x0, k, b): 
#     return L / (1 + np.exp(-k*(x-x0))) + b

# x = np.linspace(1, 40, 100)
# plt.plot(x, sigmoid(x, 2.9485, 27.507, 0.87467, 0.01047))
# plt.axvline(27.507)
# plt.show()