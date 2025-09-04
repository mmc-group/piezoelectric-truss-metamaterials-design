import torch
import numpy as np
import pandas as pd
import time
from torch.utils.data import TensorDataset, DataLoader
from torch.utils.data import random_split

import os
import sys
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from models.parameters import *
from models.model import *
from models.utils import *
from tqdm import trange
from torch.optim.lr_scheduler import StepLR

import os

torch.cuda.empty_cache()
torch.manual_seed(0)

def test(model, e_model, test_loader, saveResults, test_batch_size):

    model.eval()
    e_model.eval()

    adj_test_mse = 0.; 
    x_test_mse = 0.
    e_test_mse = 0.
    test_kld_loss = 0.

    x_test = []; adj_test = []
    x_pred = []; adj_pred = []

    with torch.no_grad():
        for adj, x, e in test_loader:
            adj = adj.to(device)
            x = x.to(device)
            e = e.to(device)

            encoded, mu, std = model.encoder(adj, x)
            if add_noise == True:
                e_input = encoded
            else:
                e_input = mu
                
            e_pred = e_model(e_input)

            adj_decoded, x_decoded = model.decoder(encoded)

            adj_mse, x_mse = recon_criterion(adj_decoded, adj), recon_criterion(x_decoded, x)
            adj_test_mse += adj_mse.item()
            x_test_mse += x_mse.item()

            test_kld = kld_loss(mu, std)
            test_kld_loss += test_kld.item()

            e_mse = stiffness_weighted_loss(e_pred, e)
            e_test_mse += e_mse.item()


    if saveResults == True:

        x_test = x.cpu().detach().numpy()
        x_pred = x_decoded.cpu().detach().numpy()

        adj_test = adj.cpu().detach().numpy()
        adj_pred = adj_decoded.cpu().detach().numpy()
        
        e_pred = e_pred.cpu().detach().numpy()
        e_test = e.cpu().detach().numpy()

        z_pred = encoded.cpu().detach().numpy()
        mu_pred = mu.cpu().detach().numpy()
        
        np.savetxt( outputFolder+'/x_test.csv', x_test, delimiter=",")
        np.savetxt( outputFolder+'/x_pred.csv', x_pred, delimiter=",")
        np.savetxt( outputFolder+'/adj_test.csv', adj_test, delimiter=",")
        np.savetxt( outputFolder+'/adj_pred.csv', adj_pred, delimiter=",")
        np.savetxt( outputFolder+'/e_test.csv', e_test, delimiter=",")
        np.savetxt( outputFolder+'/e_pred.csv', e_pred, delimiter=",")
        np.savetxt( outputFolder+'/z_pred.csv', z_pred, delimiter=",")
        np.savetxt( outputFolder+'/mu_pred.csv', mu_pred, delimiter=",")
        
    return adj_test_mse, x_test_mse, e_test_mse, test_kld

def train(epoch, KLweight):

    e_model.train()

    adj_loss_mse = 0.
    x_loss_mse = 0.
    e_loss_mse = 0.
    z_sim_loss = 0.

    data = train_prefetcher.next()
    n_batch = 0

    while data is not None:

        n_batch += 1
        if n_batch >= num_iters:
            break
        data = train_prefetcher.next()

        adj = data[0]; adj = adj.to(device)
        x = data[1]; x = x.to(device)
        e = data[2]; e = e.to(device)

        if load_pretrained_model == False:
            model.train()
            optimizer.zero_grad(set_to_none=True)
        else:
            model.eval()
            optimizer.zero_grad(set_to_none=True)

        encoded, mu, std = model.encoder(adj, x)
        if add_noise == True:
            e_input = encoded
        else:
            e_input = mu
            
        e_pred = e_model(e_input)

        adj_decoded, x_decoded = model.decoder(encoded)

        val_pred = torch.sum(adj_decoded, dim = 1)
        val = torch.sum(adj, dim = 1)

        adj_train_mse, x_train_mse = recon_criterion(adj_decoded, adj), recon_criterion(x_decoded, x)
        train_kld = kld_loss(mu, std)
        
        e_train_mse = stiffness_weighted_loss(e_pred, e)
        val_loss = recon_criterion(val, val_pred)

        if load_pretrained_model == False:
            loss = (adj_train_mse + val_loss)*a_weight + (train_kld) + x_train_mse*x_weight + e_train_mse*e_weight[epoch]
        else:
            loss = e_train_mse

        loss.backward()
        optimizer.step()
        adj_loss_mse += adj_train_mse.item()
        x_loss_mse += x_train_mse.item()
        e_loss_mse += e_train_mse.item()

    return adj_loss_mse, x_loss_mse, e_loss_mse, z_sim_loss

x_weight = 1e2
a_weight = 1
e_weight = np.ones([epochs])*1e4 # constant weight of stiffness prediction loss

model = vaeModel()
e_model = eModel()
model.to(device)
e_model.to(device)

torch.autograd.set_detect_anomaly(True)

if load_pretrained_model == True:
    model.load_state_dict(torch.load(savedModelFolder+'/best_model.pt'))
else:
    model.apply(weights_init)
e_model.apply(weights_init)

model.eval()
saveResults = None

dataset  = TensorDataset(adj_list, x, e_data)

num_train = len(dataset)
split = int(np.floor(valid_size * num_train))               # 5% of the data is used for validation
train_dataset, test_dataset = random_split(dataset = dataset, lengths = [num_train - split,split])
train_loader = DataLoader(train_dataset, batch_size = batch_size, shuffle = True)
test_loader = DataLoader(test_dataset, batch_size = test_batch_size, shuffle = True)
num_iters = len(train_loader)

if load_pretrained_model == True:
    optimizer = torch.optim.Adam(e_model.parameters(), lr = e_lr)
else:
    optimizer = torch.optim.Adam(list(model.parameters())+list(e_model.parameters()), lr = learningRate)
scheduler = StepLR(optimizer, step_size=20, gamma=0.5)

best_e_accuracy = 0.
best_adj_x_accuracy = 0.
best_accuracy_file = None
loss_file = outputFolder + '/loss.txt'
with open(loss_file, 'w') as the_file:
    the_file.write('EPOCH, Adj Train Loss, Adj Test Loss, X Train Loss, X Test Loss, e Train Loss, e Test Loss, Lr' +'\n')
the_file.close()

if best_accuracy_file is None:
    epoch_loss = [ [0 for col in range(8)] for row in range(epochs)]
    for epoch in trange(epochs):
        
        start1 = time.time()

        train_prefetcher = data_prefetcher(train_loader)
        KLweight = 1.0

        adj_loss_mse, x_loss_mse, e_loss_mse, z_sim_loss = train(epoch, KLweight)
        
        end1 = time.time()
        scheduler.step()
        
        if epoch == epochs - 1:
            saveResults = True

        adj_test_mse, x_test_mse, e_test_mse, test_kld = test(model, e_model, test_loader, saveResults, test_batch_size)
        
        epoch_loss[epoch][0] = epoch
        epoch_loss[epoch][1] = adj_test_mse
        epoch_loss[epoch][2] = x_test_mse
        epoch_loss[epoch][3] = e_test_mse

        epoch_loss[epoch][4] = adj_loss_mse
        epoch_loss[epoch][5] = x_loss_mse
        epoch_loss[epoch][6] = e_loss_mse
        epoch_loss[epoch][7] = optimizer.param_groups[0]['lr']

        if (epoch % 10 == 0.) or (epoch == epochs - 1):

            print("Epoch:", '%03d' % (epoch + 1),'/', str(epochs)\
            ,", adj train mse =", "{:.4f}".format(adj_loss_mse / len(train_loader)/batch_size)\
            ,", adj test mse =", "{:.4f}".format(adj_test_mse/len(test_loader)/test_batch_size)\
            ,", x train =", "{:.4f}".format(x_loss_mse/ len(train_loader)/batch_size )\
            ,", x test =", "{:.4f}".format(x_test_mse/len(test_loader)/test_batch_size)\
            ,", e train =", "{:.4f}".format(e_loss_mse/ len(train_loader)/batch_size )\
            ,", e test =", "{:.4f}".format(e_test_mse/len(test_loader)/test_batch_size)\
            ,", KL weight =", "{:.4f}".format(KLweight)\
            ,", time =", "{:.2f}".format(end1-start1))

            if best_e_accuracy == 0. or (e_test_mse/len(test_loader)/test_batch_size < best_e_accuracy):
                print('updating best e accuracy: previous best = {:.4f} new best = {:.4f}'.format(best_e_accuracy,
                                                                                         e_test_mse/len(test_loader)/test_batch_size))
                best_e_accuracy = e_test_mse/len(test_loader)/test_batch_size
                torch.save(e_model.state_dict(), outputFolder+'/best_e_model.pt')
                torch.save(model.state_dict(), outputFolder+'/best_model.pt')
            
model.cpu()
e_model.cpu()

model.eval()
e_model.eval()

torch.save(model.state_dict(), outputFolder+'/model.pt')
torch.save(e_model.state_dict(), outputFolder+'/e_model.pt')

loss_record = pd.DataFrame(epoch_loss, columns = ['EPOCH' , 'Adj Train Loss', 'Adj Test Loss','X Train Loss',  'X Test Loss','e Train Loss', 'e Test Loss', 'Lr'])
loss_record.to_csv(loss_file, index=False, sep=',')


