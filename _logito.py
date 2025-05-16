import string, sys, csv, os

import pandas as pd

import numpy as np

from scipy.optimize import curve_fit

import matplotlib.pyplot as plt



PY3 = sys.version_info > (3,)



class UnicodeReader:

    def __init__(self, filename, dialect=csv.excel,

                 encoding="utf-8", **kw):

        self.filename = filename

        self.dialect = dialect

        self.encoding = encoding

        self.kw = kw



    def __enter__(self):

        if PY3:

            self.f = open(self.filename, 'rt', encoding=self.encoding, newline='')

        else:

            self.f = open(self.filename, 'rb')

        self.reader = csv.reader(self.f, dialect=self.dialect, **self.kw)

        return self



    def __exit__(self, type, value, traceback):

        self.f.close()



    def next(self):

        row = next(self.reader)

        if PY3:

            return row

        return [s.decode("utf-8") for s in row]



    __next__ = next



    def __iter__(self):

        return self



class UnicodeWriter:

    def __init__(self, filename, dialect=csv.excel,

                 encoding="utf-8", **kw):

        self.filename = filename

        self.dialect = dialect

        self.encoding = encoding

        self.kw = kw



    def __enter__(self):

        if PY3:

            self.f = open(self.filename, 'wt', encoding=self.encoding, newline='')

        else:

            self.f = open(self.filename, 'wb')

        self.writer = csv.writer(self.f, dialect=self.dialect, **self.kw)

        return self



    def __exit__(self, type, value, traceback):

        self.f.close()



    def writerow(self, row):

        if not PY3:

            row = [s.encode(self.encoding) for s in row]

        self.writer.writerow(row)



    def writerows(self, rows):

        for row in rows:

            self.writerow(row)



def sigmoid(x, L ,x0, k, b):

    y = L / (1 + np.exp(-k*(x-x0)))+b

    return (y)  



def position_list():

    pos_list = []

    alpha = list(string.ascii_uppercase)[0:8]

    for i in alpha:

        for j in range(12):

            pos_list.append(i+str(j+1))

    return pos_list



def csv_write(array, file):

    with UnicodeWriter(file) as writer:

        writer.writerows(array)

    return



file_list = ['210812_R5.csv', '210813_R2.csv', '210813_R4.csv', '210813_R5.csv', '210813_R7.csv', '210813_R8.csv', '210813_R9.csv', '210813_R10.csv', '210813_R11.csv', '210814_R5.csv']



parent_dir = 'C:\\Users\\Ataturk\\Desktop\\Logistic\\2021-09-06_logito\\'

for filename in file_list:

    # filename = '210804_R5.csv'

    name_index = filename.index('.csv')

    path = os.path.join(parent_dir, filename[:-4])

    os.mkdir(path)

    print (path+'\\'+filename[:name_index]+'_logito.csv')



    dye_list = ['VIC', 'ROX', 'FAM']

    df = pd.read_csv(filename, skiprows=46)

    df = df.rename(columns={'Well Position': 'Position', 'Cycle Number': 'Cycle'})

    data_array = [['Position', 'Dye', 'Max', 'Ct-midpt', 'Growth', 'Baseline', 'Y-gain', 'PCR_eff']]





    for dye in dye_list:

        for pos in position_list():

            # print (dye, pos)

            df_new = df.copy().loc[:, ['Position', 'Cycle', dye]]

            df_new = df_new[df_new['Position']=='      '+pos]

            df_fig = df_new.loc[:, ['Cycle',dye]]

            df_fig[dye] = df_fig[dye].str.replace(',','')

            df_fig[dye] = df_fig[dye].astype(float)



            X = np.ravel(np.array(df_fig.copy().loc[:,['Cycle']]))

            y = np.ravel(np.array(df_fig.copy().loc[:,[dye]]))





            try:

                # ''' Fit to logistic function '''

                p0 = [max(y)-min(y), np.median(X), 1, min(y)] # this is an mandatory initial guess

                popt, pcov = curve_fit(sigmoid, X, y, p0, method='lm') #'dogbox')

                L, x0, k, b = popt

                Y_pred = sigmoid(X, L ,x0, k, b)

                Y_gain = np.max(Y_pred)-b

                pcr_eff = (Y_pred[3]-b)/(Y_pred[2]-b)

                data_array.append([pos, dye, L, x0, k, b, Y_gain, pcr_eff])



                ''' Figure generation '''

                # if Y_gain>100000000: #200000:

                print (pos, dye, popt, Y_gain) 

                fig = plt.figure()

                ax = fig.add_subplot(111)

                ax.scatter(X,y)

                ax.plot(X,Y_pred, color='red')

                ax.set_xlabel('Cycle number')

                ax.set_ylabel('Fluorescence')

                ax.set_title(pos + ' ' + dye, fontsize=14)

                maxpos = np.max(y)

                minpos = np.min(y)

                ax.text(-5, 1.00*(maxpos-minpos)+minpos, 'Y_gain: '+str("{:.3e}".format(Y_gain)), fontsize=14)

                ax.text(-5, 0.95*(maxpos-minpos)+minpos, 'PCR eff: '+str("{:.2f}".format(pcr_eff)), fontsize=14)

                ax.text(-5, 0.90*(maxpos-minpos)+minpos, 'Ct at midpt: '+str("{:.2f}".format(popt[1])), fontsize=14)

                ax.text(-5, 0.85*(maxpos-minpos)+minpos, 'Growth: '+str("{:.3e}".format(popt[0])), fontsize=14)

                ax.text(-5, 0.80*(maxpos-minpos)+minpos, 'Baseline: '+str("{:.3e}".format(popt[3])), fontsize=14)

                plt.savefig(path+'\\'+filename[:name_index]+dye+'_'+pos, bbox_inches='tight')

                plt.close()



            except Exception as e:

                print('error', pos, dye, e)

                pass



    # ''' save logistic function parameters '''

    csv_write(data_array, path+'\\'+filename[:name_index]+'_logito.csv')
