from fidelity_function import *
from qutip import *
import numpy as np
from matplotlib import pyplot as plt

kwargs = {'I1_p': 0, 'Q1_p':50 , 'I2_p': 0, 'Q2_p': 0 }
fid_paper = fidelity_fn_X(**kwargs)