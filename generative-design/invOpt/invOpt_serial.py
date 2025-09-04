import torch
import numpy as np
import time
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

import os
import sys
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from models.parameters import *
from models.model import *
from models.utils import *


from tqdm import tqdm
import networkx as nx
from networkx import connected_components
import os.path

torch.autograd.set_detect_anomaly(True)

torch.manual_seed(0)
torch.backends.cudnn.benchmark = True


model = vaeModel()
e_model = eModel()

model.load_state_dict(torch.load('../results/best_model.pt', map_location=torch.device(device)))
e_model.load_state_dict(torch.load('../results/best_e_model.pt', map_location=torch.device(device)))

model.eval()
e_model.eval()

model.to(device)
e_model.to(device)
ptb_mask = pd.read_csv(folder+'/ptb_mask.csv', delimiter = ",", header = None).to_numpy()

a_row, a_col = np.triu_indices(numNodes)
x_row, x_col = np.nonzero(ptb_mask)

relu_func = nn.ReLU()

# Get task
opt_target_code = sys.argv[1]
seed = sys.argv[2]
np.random.seed(int(seed))
torch.manual_seed(int(seed))
w = 10.

def recon_criterion(pred, opt_target_code):
    pred = unnormalize(pred)
    if opt_target_code == '+e31-e33':
        return -pred[:,2]+torch.abs(pred[:,4])
    elif opt_target_code == 'custom':
        # Define custom objective
        return
    else:
        print('Missing / Unknown objective')
    
outputFolder = '../results'
inverseSaveFolder = outputFolder+'/opt/'
inverseSaveFolder += opt_target_code+'_seed='+seed

os.makedirs(inverseSaveFolder,exist_ok=True)

dataset  = TensorDataset(adj_list, x, e_data)
num_train = len(dataset)
data_loader = DataLoader(dataset, batch_size = num_train, shuffle = False)

optimization_method = ['Adam']
opt_epoch = 30000
patience = 25000
num_sample = 1
Adam_lr = 1e-4
LBFGS_lr = 1e-2
trace_back = True
num_workers = 1

start = time.time()
var_dict = {}
torch.manual_seed(int(seed))
np.random.seed(int(seed))

if os.path.exists(inverseSaveFolder + "/train_std.csv"):
    train_z = pd.read_csv(inverseSaveFolder+"/train_z.csv", delimiter = ",", header = None).to_numpy()
    train_mu = pd.read_csv(inverseSaveFolder+"/train_mu.csv", delimiter = ",", header = None).to_numpy()
    train_std = pd.read_csv(inverseSaveFolder+"/train_std.csv", delimiter = ",", header = None).to_numpy()
else:
    train_z = []
    train_mu = []
    train_std = []
    for adj, x, c in data_loader:
        adj = adj.to(device)
        x = x.to(device)
        c = c.to(device)
        encoded, mu, std = model.encoder(adj, x)
        train_z.extend(encoded.detach().cpu().numpy())
        train_mu.extend(mu.detach().cpu().numpy())
        train_std.extend(std.detach().cpu().numpy())

dtype = torch.float
train_z = np.array(train_z)

print("**************************************")
print("***** Inverse design target ", opt_target_code," *****")
print("**************************************") 

try:
    os.makedirs(inverseSaveFolder+'/',exist_ok=True)
except OSError:
    print ("Creation of the directory failed")
else:
    print ('Successfully created the directory ' + inverseSaveFolder)

stiffness_vec = e_data.numpy()
moduli = pd.read_csv(folder+'/e_data.csv', delimiter = ",", header = None).to_numpy()
moduli = unnormalize(torch.tensor(moduli)).numpy()
loss = 0.

if opt_target_code == '+e31-e33':
    best_val = np.max(moduli[:,2]-np.abs(moduli[:,4]))
    prop_arr = -moduli[:,2]+np.abs(moduli[:,4])
    loss += prop_arr
elif opt_target_code == 'custom':
    # Define custom objective
    pass
else:
    print('Missing / Unknown objective')

moduli = normalize(torch.tensor(moduli)).numpy()

print("Best objective value in dataset: ",opt_target_code," --- ",best_val)
    
initial_idx = np.argsort(loss)[:int(num_sample)]
print("Number of initial guesses = ", len(initial_idx))
initial_adj = adj_list[initial_idx,:].detach().clone().float().to(device)
initial_x = x[initial_idx,:].detach().clone().float().to(device)

initial_z, initial_mu, initial_std = model.encoder(initial_adj, initial_x)

var_dict = {}
updated_z = np.zeros(initial_z.shape)
predicted_e= np.zeros([len(initial_z), paramDim])
optimization = optimization_method[0]

def main(opt_target_code):

    best_guess_acc = 100.
    best_acc = 100.

    z_trace_np = np.zeros((num_sample, opt_epoch, latent_dim))
    e_trace_np = np.zeros((num_sample, opt_epoch, moduli.shape[1]))
    x_trace_np = np.zeros((num_sample, opt_epoch, initial_x.shape[1]))
    adj_trace_np = np.zeros((num_sample, opt_epoch, initial_adj.shape[1]))
    updated_z_np = np.zeros((num_sample, latent_dim))
    predicted_e_np = np.zeros((num_sample, moduli.shape[1]))
    
    
    for idx in range(len(initial_z)):
        es=0
        initial_guess_z = (initial_z[idx,:]).reshape([1, initial_z.shape[-1]])
        initial_guess_z = torch.tensor(initial_guess_z.clone(), device=device).float().requires_grad_()

        if optimization == 'Adam':
            opt = torch.optim.Adam([initial_guess_z], lr = Adam_lr)

        elif optimization == 'NAdam':
            opt = torch.optim.NAdam([initial_guess_z], lr = Adam_lr)
        
        elif optimization == 'LBFGS':
            opt = torch.optim.LBFGS([initial_guess_z],lr = LBFGS_lr)

        progress_bar = tqdm(total=opt_epoch, desc='Optimization')

        for e in range(opt_epoch):

            inv_adj_recon, inv_x_recon = model.decoder(initial_guess_z)
            inv_x_recon = relu_func(inv_x_recon - 5e-3)
            inv_adj_recon = relu_func(inv_adj_recon - 0.1)
            inv_z_recon, inv_mu_recon, _ = model.encoder(inv_adj_recon, inv_x_recon)
          
            inv_e_pred = e_model(inv_mu_recon)
            loss = 0.
            loss += recon_criterion(inv_e_pred, opt_target_code=opt_target_code) 

            loss.backward(retain_graph = True)
            opt.step(lambda: recon_criterion(inv_e_pred, opt_target_code=opt_target_code)) 
            valid_flag = validity_check(inv_adj_recon,inv_x_recon)
            if valid_flag:
                z_trace_np[idx,e,:] = inv_z_recon.cpu().detach().numpy()
                e_trace_np[idx,e,:] = inv_e_pred.cpu().detach().numpy()
                x_trace_np[idx,e,:] = inv_x_recon.cpu().detach().numpy()
                adj_trace_np[idx,e,:] = inv_adj_recon.cpu().detach().numpy()

            if loss.item() < best_acc:

                if valid_flag:
                    best_acc = loss.item()
                    best_z = inv_z_recon.cpu().detach().numpy()
                    best_e = inv_e_pred.cpu().detach().numpy()
                    best_x = inv_x_recon.cpu().detach().numpy()
                    best_adj = inv_adj_recon.cpu().detach().numpy()
                    es = 0
                    
            else:
                es += 1
                if es > patience:
                    break
            
            progress_bar.set_description(f'Objective : {best_acc:.4f}')
            progress_bar.update(1)

        try:
            updated_z_np[idx,:] = best_z
            predicted_e_np[idx,:] = best_e
        except:
            updated_z_np[idx,:] = inv_z_recon.cpu().detach()
            predicted_e_np[idx,:] = inv_e_pred.cpu().detach()
            continue

        if best_acc < best_guess_acc:
            best_inv_pred = idx
            best_guess_acc = best_acc

    print("Predicted best id: ", best_inv_pred)
    os.makedirs(inverseSaveFolder,exist_ok=True)
    np.savetxt(inverseSaveFolder+'/best_e.txt', predicted_e_np[best_inv_pred,:].flatten())
    np.savetxt(inverseSaveFolder+'/best_z.csv', updated_z_np[best_inv_pred,:], delimiter = ",")
    np.save(inverseSaveFolder+'/best_x', best_x)
    np.save(inverseSaveFolder+'/best_adj', best_adj)
    e_unnorm = unnormalize(torch.tensor(predicted_e_np[best_inv_pred:best_inv_pred+1,:]))
    np.savetxt(inverseSaveFolder+'/best_e_unnorm.txt', e_unnorm.flatten())
    np.savetxt(inverseSaveFolder+'/best_e_dataset.txt', moduli[initial_idx[0]].flatten())
    np.save(inverseSaveFolder+'/z_trace_best', z_trace_np[best_inv_pred,:,:])
    np.save(inverseSaveFolder+'/x_trace_best', x_trace_np[best_inv_pred,:,:])
    np.save(inverseSaveFolder+'/adj_trace_best', adj_trace_np[best_inv_pred,:,:])
    np.save(inverseSaveFolder+'/e_trace_best', unnormalize(torch.tensor(e_trace_np[best_inv_pred,:,:])).numpy())


def validity_check(inv_adj_recon,inv_x_recon):
    inv_adj_recon = inv_adj_recon.detach().cpu().numpy()
    inv_x_recon = inv_x_recon.detach().cpu().numpy()

    a_recon = adj_vec2array(inv_adj_recon, a_row, a_col)

    new_base_lattice = np.array(a_recon)
    ex = new_base_lattice
    ex[np.triu_indices(numNodes)] = 0.
    row = np.nonzero(ex)[0]
    col = np.nonzero(ex)[1]
    g = nx.Graph()
    for i in range(row.shape[0]):
        g.add_edge(row[i], col[i])
    num_connected = len(list(list(connected_components(g))))
    if num_connected == 0 or num_connected > 1:
        return False
    elif num_connected == 1:
        return True
    else:
         print('Error')
        
    
import pdb
import sys
import traceback

if __name__ == '__main__':
    try:
        main(opt_target_code)
    except:
        extype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)

