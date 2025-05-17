import os
import matplotlib.pyplot as plt
from Ripper import Ripper

fpath = os.path.join("apldbio", "sds", "analysis_result.txt")
ripper = Ripper(fpath)

with open(fpath) as fhandle:
    data = ripper._default_parse_result_file(fhandle)
    key = "A9"
    x0 = ripper.fit_logito(data[key]["Delta Rn values"])
    # plt.plot(data[key]["Delta Rn values"])
    # plt.axvline(x0)
    # plt.show()