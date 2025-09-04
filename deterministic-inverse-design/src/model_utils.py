import torch
import torch.nn.functional as F
from parameters import *
from src.voigt_rotation import *
import pickle, io

# unpickle object also with a CPU-only machine, see issue: https://github.com/pytorch/pytorch/issues/16797
class CPU_Unpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == 'torch.storage' and name == '_load_from_bytes':
            return lambda b: torch.load(io.BytesIO(b), map_location='cpu')
        else: return super().find_class(module, name)

def getActivation(activ):
    if(activ == 'relu'):
        sigma = torch.nn.ReLU()
    elif(activ == 'tanh'):
        sigma = torch.nn.Tanh()
    elif(activ == 'sigmoid'):
        sigma = torch.nn.Sigmoid()
    elif(activ == 'leaky'):
        sigma = torch.nn.LeakyReLU()
    elif(activ == 'softplus'):
        sigma = torch.nn.Softplus()
    elif(activ == 'logsigmoid'):
        sigma = torch.nn.LogSigmoid()
    elif(activ == 'elu'):
        sigma = torch.nn.ELU()
    elif(activ == 'gelu'):
        sigma = torch.nn.GELU()
    elif(activ == 'none'):
        sigma = torch.nn.Identity()
    else:
        raise ValueError('Incorrect activation function')
    return sigma

def createNN(inputDim,arch,outputDim,bias=True):
    model = torch.nn.Sequential()
    currDim = inputDim
    layerCount = 1
    activCount = 1
    for i in range(len(arch)):
        if(type(arch[i]) == int):
            model.add_module('layer '+str(layerCount),torch.nn.Linear(currDim,arch[i],bias=bias))
            currDim = arch[i]
            layerCount += 1
        elif(type(arch[i]) == str):
            model.add_module('activ '+str(activCount),getActivation(arch[i]))
            activCount += 1
    model.add_module('layer '+str(layerCount),torch.nn.Linear(currDim,outputDim,bias=bias))
    return model

def softmax(input, t):
    return F.log_softmax(input/t, dim=1)

def gumbel(input, t):
    return F.gumbel_softmax(input, tau=t, hard=True, eps=1e-10, dim=1)

def assemble_F2_features(C_ort,e_rotated,R1,V,C_ort_scaling,e_ort_scaling,C_scaling,e_scaling,method=None):
    # scale C_ort to its original range
    C_ort_unscaled = C_ort_scaling.unnormalize(C_ort)
    # rotate C_ort (directly in Voigt notation)
    C_tilde = direct_rotate(C_ort_unscaled,R1,orthotropic=True,method=method)
    C_tilde = C_scaling.normalize(C_tilde)
    return torch.cat((C_tilde,e_rotated,V),dim=1)

def invModel_output(G1,G2,input,t,activation):
    # continuous params: [stretch1, stretch2, stretch3, rot_stretch1, rot_stretch2, rot_stretch3, theta, rot_ax1, rot_ax2]
    topology1,topology2,topology3,rep1,rep2,rep3 = torch.split(G1(input), [7,7,7,2,2,2], dim=1)
    m = getActivation('sigmoid')
    if(activation == 'one-hot'):
        # enforce one-hot encoding by small temperature
        t = 1.e-6
    if(activation == 'softmax' or activation == 'one-hot'):
        topology = torch.cat((softmax(topology1,t),softmax(topology2,t),softmax(topology3,t),softmax(rep1,t),softmax(rep2,t),softmax(rep3,t)), dim=1)
    elif(activation == 'gumbel'):
        topology1,topology2,topology3,rep1,rep2,rep3 = softmax(topology1,t),softmax(topology2,t),softmax(topology3,t),softmax(rep1,t),softmax(rep2,t),softmax(rep3,t)
        topology = torch.cat((gumbel(topology1,t),gumbel(topology2,t),gumbel(topology3,t),gumbel(rep1,t),gumbel(rep2,t),gumbel(rep3,t)), dim=1)
    else:
        raise ValueError('Incorrect activation function')

    features = torch.cat((topology, input), dim=1)
    rho_U, V, rot1, rot2 = torch.split(G2(features), [4,3,6,6], dim=1)
    # scale to [0,1] using sigmoid
    rho_U, V = m(rho_U), m(V)

    return rho_U, V, rot1, rot2, topology

def rotate_C(C_in,R,C_in_scaling,C_out_scaling,method=None):
    temp = C_in_scaling.unnormalize(C_in)
    temp = direct_rotate(temp,R,method=method)
    C = C_out_scaling.normalize(temp)
    
    return C


def rotate_e(e_in,R,e_in_scaling,e_out_scaling,method=None):

    temp = e_in_scaling.unnormalize(e_in)
    temp = direct_rotate_e(temp,R)
    e = e_out_scaling.normalize(temp)
    return e

def assemble_C_matrix(C):

    C_matrix = torch.zeros((len(C),6,6))

    ## Assemble C_matrix
    C_matrix[:,0,0] = C[:,0]
    C_matrix[:,0,1] = C[:,1]
    C_matrix[:,0,2] = C[:,2]
    C_matrix[:,0,3] = C[:,3]
    C_matrix[:,0,4] = C[:,4]
    C_matrix[:,0,5] = C[:,5]

    C_matrix[:,1,0] = C[:,1]
    C_matrix[:,1,1] = C[:,6]
    C_matrix[:,1,2] = C[:,7]
    C_matrix[:,1,3] = C[:,8]
    C_matrix[:,1,4] = C[:,9]
    C_matrix[:,1,5] = C[:,10]

    C_matrix[:,2,0] = C[:,2]
    C_matrix[:,2,1] = C[:,7]
    C_matrix[:,2,2] = C[:,11]
    C_matrix[:,2,3] = C[:,12]
    C_matrix[:,2,4] = C[:,13]
    C_matrix[:,2,5] = C[:,14]

    C_matrix[:,3,0] = C[:,3]
    C_matrix[:,3,1] = C[:,8]
    C_matrix[:,3,2] = C[:,12]
    C_matrix[:,3,3] = C[:,15]
    C_matrix[:,3,4] = C[:,16]
    C_matrix[:,3,5] = C[:,17]

    C_matrix[:,4,0] = C[:,4]
    C_matrix[:,4,1] = C[:,9]
    C_matrix[:,4,2] = C[:,13]
    C_matrix[:,4,3] = C[:,16]
    C_matrix[:,4,4] = C[:,18]
    C_matrix[:,4,5] = C[:,19]

    C_matrix[:,5,0] = C[:,5]
    C_matrix[:,5,1] = C[:,10]
    C_matrix[:,5,2] = C[:,14]
    C_matrix[:,5,3] = C[:,17]
    C_matrix[:,5,4] = C[:,19]
    C_matrix[:,5,5] = C[:,20]

    return C_matrix

def assemble_e_matrix(e):

    e_matrix = torch.zeros((len(e),3,6))

    ## Assemble e_matrix
    e_matrix[:,0,0] = e[:,0]
    e_matrix[:,0,1] = e[:,1]
    e_matrix[:,0,2] = e[:,2]
    e_matrix[:,0,3] = e[:,3]
    e_matrix[:,0,4] = e[:,4]
    e_matrix[:,0,5] = e[:,5]

    e_matrix[:,1,0] = e[:,6]
    e_matrix[:,1,1] = e[:,7]
    e_matrix[:,1,2] = e[:,8]
    e_matrix[:,1,3] = e[:,9]
    e_matrix[:,1,4] = e[:,10]
    e_matrix[:,1,5] = e[:,11]

    e_matrix[:,2,0] = e[:,12]
    e_matrix[:,2,1] = e[:,13]
    e_matrix[:,2,2] = e[:,14]
    e_matrix[:,2,3] = e[:,15]
    e_matrix[:,2,4] = e[:,16]
    e_matrix[:,2,5] = e[:,17]

    return e_matrix

def compute_d_coeff(C_matrix,e_matrix,full_matrix=False):

    S = torch.linalg.inv(C_matrix)
    d_matrix = torch.matmul(e_matrix,S)

    if full_matrix:

        d = torch.zeros((len(d_matrix),18))

        d[:,0] = d_matrix[:,0,0]
        d[:,1] = d_matrix[:,0,1]
        d[:,2] = d_matrix[:,0,2]
        d[:,3] = d_matrix[:,0,3]
        d[:,4] = d_matrix[:,0,4]
        d[:,5] = d_matrix[:,0,5]

        d[:,6] = d_matrix[:,1,0]
        d[:,7] = d_matrix[:,1,1]
        d[:,8] = d_matrix[:,1,2]
        d[:,9] = d_matrix[:,1,3]
        d[:,10] = d_matrix[:,1,4]
        d[:,11] = d_matrix[:,1,5]

        d[:,12] = d_matrix[:,2,0]
        d[:,13] = d_matrix[:,2,1]
        d[:,14] = d_matrix[:,2,2]
        d[:,15] = d_matrix[:,2,3]
        d[:,16] = d_matrix[:,2,4]
        d[:,17] = d_matrix[:,2,5]

    else:
        d = torch.zeros((len(d_matrix),1))
        d[:,0] = d_matrix[:,2,2]
        # d[:,1] = d_matrix[:,1,3]
        # d[:,2] = d_matrix[:,2,2]

    return d
