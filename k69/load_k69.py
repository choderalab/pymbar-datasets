import numpy
 
def load_dataset():
    npzfile = numpy.load('k69_dataset.npz')
    #u_kln matrix
    u_kln = npzfile['u']
    #N_k uncorrelated samples per state
    N_k = npzfile['N_k']
    #Converged mbar weights to check against
    f_k = npzfile['f_k']
    return u_kln, N_k, f_k


u_kln, N_k, f_k = load_dataset()
import pymbar

mbar = pymbar.MBAR(u_kln, N_k)
pymbar.testsystems.pymbar_datasets.save("k69", mbar.u_kn, mbar.N_k, least_significant_digit=8)
