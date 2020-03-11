import numpy as np
import math
# Credits to Jens Glaser
#Implements the method of H. Flyvbjerg and H. G. Petersen
#  doi: 10.1063/1.457480

class BlockAverage:
    def __init__(self, ydata, minlen=16):
        block = np.array(ydata)
        self.err_est = []
        self.err_err = []
        self.n = []
        self.num_samples = []
        i = 0
        while (len(block) >= int(minlen)):
            e = 1.0/(len(block)-1)*np.var(block)
            self.err_est.append(math.sqrt(e))
            self.err_err.append(math.sqrt(e/2.0/(len(block)-1)))
            self.n.append(i)
            self.num_samples.append(len(block))
            block_l = block[1:]
            block_r = block[:-1]
            block = 1.0/2.0*(block_l + block_r)
            block = block[::2]
            i += 1

        # convert to np arrays
        self.err_est = np.array(self.err_est)
        self.err_err = np.array(self.err_err)
        self.num_samples = np.array(self.num_samples)

    def get_hierarchical_errors(self):
        return (self.n, self.num_samples, self.err_est, self.err_err)

    def get_error_estimate(self,relsigma=1.0):
        i = self.n[-1]
        while True:
            # weighted error average
            avg_err = np.sum(self.err_est[i:]/self.err_err[i:]/self.err_err[i:])
            avg_err /= np.sum(1.0/self.err_err[i:]/self.err_err[i:])

            sigma = self.err_err[i]
            cur_err = self.err_est[i]
            delta = abs(cur_err - avg_err)
            if (delta > relsigma*sigma or i == 0):
                i += 1
                break
            i -= 1

        # compute average error in plateau region
        avg_err = np.sum(self.err_est[i:]/self.err_err[i:]/self.err_err[i:])
        avg_err /= np.sum(1.0/self.err_err[i:]/self.err_err[i:])
        return (i,avg_err)
