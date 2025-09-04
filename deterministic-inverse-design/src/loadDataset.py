import torch
from torch.utils.data import TensorDataset
import pickle
import numpy as np
import pandas as pd
from parameters import *
from src.normalization import Normalization
from src.voigt_rotation import *
from src.model_utils import CPU_Unpickler, assemble_C_matrix, assemble_e_matrix, compute_d_coeff

def exportTensor(name,data,cols, header=True):
    df=pd.DataFrame.from_records(data.detach().numpy())
    if(header):
        df.columns = cols
    print(name)
    df.to_csv(name+".csv", header=header, index=False)

def exportList(name,data):
    arr=np.array(data)
    np.savetxt(name+".csv", [arr], delimiter=',')

def getNormalization(save_normalization=False):

    data_inp = pd.read_csv(dataPath_input)
    data_out = pd.read_csv(dataPath_output)

    # check for NaNs
    assert not data_inp.isnull().values.any()
    assert not data_out.isnull().values.any()

    F1_features = torch.tensor(data_inp[features_names].values)
    R2 = torch.tensor(data_inp[R2_names].values)
    V = torch.tensor(data_inp[V_names].values)
    e = torch.tensor(data_out[e_names].values)

    F1_features_scaling = Normalization(F1_features,features_types,features_scaling_strategy)
    V_scaling = Normalization(V,V_types,V_scaling_strategy)
    e_scaling = Normalization(e,e_types,e_scaling_strategy)


    # should only be activated if framework is retrained with different dataset
    if save_normalization:
        with open('src/normalization/F1_features_scaling.pickle', 'wb') as file_:
            pickle.dump(F1_features_scaling, file_, -1)
        with open('src/normalization/V_scaling.pickle', 'wb') as file_:
            pickle.dump(V_scaling, file_, -1)
        with open('src/normalization/e_scaling.pickle', 'wb') as file_:
            pickle.dump(e_scaling, file_, -1)

    return F1_features_scaling, e_scaling, V_scaling

def getSavedNormalization():

    F1_features_scaling = CPU_Unpickler(open("src/normalization/F1_features_scaling.pickle", "rb", -1)).load()
    V_scaling = CPU_Unpickler(open("src/normalization/V_scaling.pickle", "rb", -1)).load()
    e_scaling = CPU_Unpickler(open("src/normalization/e_scaling.pickle", "rb", -1)).load()

    return  F1_features_scaling, e_scaling, V_scaling

def getDataset(F1_features_scaling, e_scaling, V_scaling, return_individual=False, normalize=True):

    data_inp = pd.read_csv(dataPath_input)
    data_out = pd.read_csv(dataPath_output)

    print('Input data:                ',data_inp.shape)
    print('Output data:               ',data_out.shape)

    # check for NaNs
    assert not data_inp.isnull().values.any()
    assert not data_out.isnull().values.any()

    F1_features = torch.tensor(data_inp[features_names].values)
    R1 = torch.tensor(data_inp[R1_names].values)
    R2 = torch.tensor(data_inp[R2_names].values)
    V = torch.tensor(data_inp[V_names].values)
    e = torch.tensor(data_out[e_names].values)
    
    if normalize == True:
        F1_features = F1_features_scaling.normalize(F1_features)
        V = V_scaling.normalize(V)
        e = e_scaling.normalize(e)

    if return_individual == True:
        return F1_features.float(), R1.float(), V.float(),\
            R2.float(), e.float()

    dataset =  TensorDataset(F1_features.float(), R1.float(), V.float(), R2.float(), e.float())
    l1 = round(len(dataset)*traintest_split)
    l2 = len(dataset) - l1
    print('train/test: ',[l1,l2],'\n\n')

    train_set, test_set = torch.utils.data.random_split(dataset, [l1,l2], generator=torch.Generator().manual_seed(42))
    return train_set, test_set