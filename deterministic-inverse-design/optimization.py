#------------------------------------------------------Load libraries-------------------------------------------------------
import sys
sys.path.insert(0, "../")
import numpy as np
import os
import torch
import pandas as pd
from torch import optim
import torch
from parameters import *
from src.loadDataset import *
from src.model_utils import *
import shutil
print(device)


objective = sys.argv[1] # ['max-fom','max-e31-min-e33',]
if objective is None:
    print("Provide objective argument: ['max-fom','max-e15','max-all'] or define custom objective in source code.")
else:
    print("Optimization objective: ", objective)

#---------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------Define functions for optimization----------------------------------------------
def convert_lattice_types(lattice_types):
    """ Decodes predicted softmax lattice types into array of integers."""
    lattice_types_int = torch.trunc(7*lattice_types)
    return lattice_types_int

def convert_lattice_reps(lattice_reps):
    """ Decodes predicted softmax lattice repetitions into array of integers."""
    lattice_reps_int = torch.trunc(2*lattice_reps)
    return lattice_reps_int

def get_features(solution: torch.Tensor, epoch: int):
    """Extract features from ML prediction at each epoch

    Args:
        solution (torch.Tensor): The ML prediction at current epoch of the optimization loop.
        epoch (int): The current epoch of the optimization loop

    Returns:
        F1_features_pred (torch.Tensor): The predicted topological features of the unit cell.
        R1_pred_6D: (torch.Tensor): The first rotation tensor in 6D representation.
        V_pred (torch.Tensor): The eigenvalues of the second stretch tensor.
        R2_pred_6D (torch.Tensor): The second rotation tensor in 6D representation.
    """
    t=100. # Gumbel-Softmax temeperature
    if epoch == 0:
        # F1 features
        F1_features_pred = solution[0:31]
        # R1 & R2
        R1_pred_6D = solution[31:37]
        R2_pred_6D = solution[31:37]
        # V2
        V_pred = solution[43:46]

    else:
        # F1 features
        F1_features_pred = solution[0:31]
        relative_density = torch.sigmoid(F1_features_pred[0:1])
        if limit_stretch == True:
            U = 0.5*torch.sigmoid(F1_features_pred[1:4])
        else:
            U = torch.sigmoid(F1_features_pred[1:4])

        # Convert to one hot
        lattice_types_tmp_1 = torch.nn.functional.gumbel_softmax(F1_features_pred[4:11], dim=0, tau=t, hard=True, eps=1e-10).view(1,-1)
        lattice_types_tmp_2 = torch.nn.functional.gumbel_softmax(F1_features_pred[11:18], dim=0, tau=t, hard=True, eps=1e-10).view(1,-1)
        lattice_types_tmp_3 = torch.nn.functional.gumbel_softmax(F1_features_pred[18:25], dim=0, tau=t, hard=True, eps=1e-10).view(1,-1)
        lattice_types_one_hot = torch.cat((lattice_types_tmp_1,lattice_types_tmp_2,lattice_types_tmp_3),dim=1)

        lattice_reps_tmp_1 = torch.nn.functional.gumbel_softmax(F1_features_pred[25:27], dim=0, tau=t, hard=True, eps=1e-10).view(1,-1)
        lattice_reps_tmp_2 = torch.nn.functional.gumbel_softmax(F1_features_pred[27:29], dim=0, tau=t, hard=True, eps=1e-10).view(1,-1)
        lattice_reps_tmp_3 = torch.nn.functional.gumbel_softmax(F1_features_pred[29:31], dim=0, tau=t, hard=True, eps=1e-10).view(1,-1)
        lattice_reps_one_hot = torch.cat((lattice_reps_tmp_1,lattice_reps_tmp_2,lattice_reps_tmp_3),dim=1)

        F1_features_pred = torch.cat((relative_density.view(1,-1),U.view(1,-1),lattice_types_one_hot.view(1,-1),lattice_reps_one_hot.view(1,-1)),dim=1)
        # R1 & R2
        R1_pred_6D = torch.tanh(solution[31:37])
        R2_pred_6D = torch.tanh(solution[37:43])
        # V2
        if limit_stretch == True:
            V_pred = 0.5*torch.sigmoid(solution[43:46]).view(1,-1)
        else:
            V_pred = torch.sigmoid(solution[43:46]).view(1,-1)

    return F1_features_pred,R1_pred_6D,V_pred,R2_pred_6D

def loss_func(e_matrix: torch.Tensor, objective: str) -> torch.Tensor:

    if objective == 'max-fom':
        
        loss = -(e_matrix[0,2,0] + e_matrix[0,2,1] + e_matrix[0,2,2]) \
                 + magnitude_scaler*((e_matrix[0,2,0]-e_matrix[0,2,1])**2  \
                 + (e_matrix[0,2,0]-e_matrix[0,2,2])**2 + (e_matrix[0,2,1]-e_matrix[0,2,2])**2)
        
        loss = loss - shrinkage_scaler*(-e_matrix[:,0,0]**2 - e_matrix[:,0,1]**2 - e_matrix[:,0,2]**2 \
                                     -e_matrix[:,0,3]**2 - e_matrix[:,0,5]**2 - e_matrix[:,1,0]**2 \
                                     -e_matrix[:,1,1]**2 - e_matrix[:,1,2]**2 - e_matrix[:,1,4]**2 \
                                     -e_matrix[:,1,5]**2 - e_matrix[:,2,3]**2 - e_matrix[:,2,4]**2 \
                                     -e_matrix[:,2,5]**2)

    elif objective == 'e-15':
        loss = -e_matrix[0,0,4]
        loss = loss - shrinkage_scaler*(-e_matrix[:,0,0]**2 - e_matrix[:,0,1]**2 - e_matrix[:,0,2]**2 \
                                        -e_matrix[:,0,3]**2 - e_matrix[:,0,5]**2 - e_matrix[:,1,0]**2 \
                                        -e_matrix[:,1,1]**2 - e_matrix[:,1,2]**2 - e_matrix[:,1,3]**2 \
                                        -e_matrix[:,1,4]**2 - e_matrix[:,1,5]**2 - e_matrix[:,2,1]**2 \
                                        -e_matrix[:,2,2]**2 - e_matrix[:,2,3]**2 - e_matrix[:,2,4]**2 \
                                        - e_matrix[:,2,5]**2)

    elif objective == 'max-all':
        loss = -torch.sum(e_matrix)

    elif objective == 'custom':
        # Create custom objective here.
        return
    
    return loss


def get_loss(model: torch.nn.Module, e_scaling: Normalization, predicted_parameters: torch.Tensor, epoch: int):
    """ Obtain the loss value for the current epoch.

    Args:
        predicted_parameters (torch.Tensor): The ML prediction at current epoch of the optimization loop.
        epoch (int): The current epoch of the optimization loop

    Returns:
        loss (torch.Tensor): Current loss
    """

    F1_features_pred,R1_pred_6D,V_pred,R2_pred_6D = get_features(predicted_parameters, epoch)

    e_pred = model(torch.cat((F1_features_pred.view(1,-1),R1_pred_6D.view(1,-1),V_pred.view(1,-1),R2_pred_6D.view(1,-1)),dim=1))
    e_pred = e_scaling.unnormalize(e_pred)
    e_pred_matrix = assemble_e_matrix(e_pred.view(1,-1))

    loss = loss_func(e_pred_matrix, objective)
    
    return loss

def delete_previous_results(objective: str):
    """Delete previous results based on the task type."""
    
    file_path = 'results/'+objective+'/'
    if os.path.isdir(file_path):
        shutil.rmtree(file_path)

    os.makedirs(file_path)

def save_trace(predicted_parameters: torch.Tensor, F1_features_scaling: Normalization, V_scaling: Normalization, epoch: int, num_guess: int):
    """ Writing output files at each epoch to save the trace of the optimization.

    Args:
        predicted_parameters (torch.Tensor): The predicted parameters at the current epoch.
        F1_features_scaling (Normalization): The Normalization instance for the predicted topological parameters.
        V_scaling (Normalization): The Normalization instance for the second stretch tensor.
        epoch (int): The current epoch.
        num_guess (int): The index of the current initial guess.
    """

    F1_features_pred,R1_pred_6D,V_pred,R2_pred_6D = get_features(predicted_parameters, epoch=1)

    relative_density = torch.mul(F1_features_pred[0][0:1], F1_features_scaling.max[0] - F1_features_scaling.min[0]) + F1_features_scaling.min[0]
    U1 = torch.mul(F1_features_pred[0][1:2], F1_features_scaling.max[1]-F1_features_scaling.min[1]) + F1_features_scaling.min[1]
    U2 = torch.mul(F1_features_pred[0][2:3], F1_features_scaling.max[2]-F1_features_scaling.min[2]) + F1_features_scaling.min[2]
    U3 = torch.mul(F1_features_pred[0][3:4], F1_features_scaling.max[3]-F1_features_scaling.min[3]) + F1_features_scaling.min[3]
    U = torch.cat((U1.view(1,-1),U2.view(1,-1),U3.view(1,-1)),dim=1)
    
    lattice_types_tmp = F1_features_pred[0][4:25]
    lattice_types = torch.zeros(3)

    for i in range(3):
        lattice_types[i] = torch.argmax(lattice_types_tmp[i*7:(i+1)*7])
    
    lattice_reps_tmp = F1_features_pred[0][25:31]
    lattice_reps = torch.zeros(3)

    for i in range(3):
        lattice_reps[i] = torch.argmax(lattice_reps_tmp[i*2:(i+1)*2])

    R1_3d = rot6DToAngleAxis(R1_pred_6D.view(1,-1))
    R2_3d = rot6DToAngleAxis(R2_pred_6D.view(1,-1))

    V1 = torch.mul(V_pred[0:1], V_scaling.max[0]-V_scaling.min[0]) + V_scaling.min[0]
    V2 = torch.mul(V_pred[1:2], V_scaling.max[1]-V_scaling.min[1]) + V_scaling.min[1]
    V3 = torch.mul(V_pred[2:3], V_scaling.max[2]-V_scaling.min[2]) + V_scaling.min[2]
    V = torch.cat((V1.view(1,-1),V2.view(1,-1),V3.view(1,-1)),dim=1)

    header=['relative_density','U1','U2','U3','lattice_type1','lattice_type2','lattice_type3','lattice_rep1','lattice_rep2','lattice_rep3',\
            'R1_theta','R1_rot_ax1','R1_rot_ax2','R2_theta','R2_rot_ax1','R2_rot_ax2','V1','V2','V3']
    
    input_parameter_trace = torch.cat((relative_density.view(1,-1).to(device),U.view(1,-1).to(device),lattice_types.view(1,-1).to(device),lattice_reps.view(1,-1).to(device),R1_3d.view(1,-1).to(device),\
                            R2_3d.view(1,-1).to(device),V.view(1,-1).to(device)),dim=1)
    # breakpoint()
    input_parameter_trace = pd.DataFrame(input_parameter_trace.detach().cpu().numpy())
    if epoch == 0:
        input_parameter_trace.to_csv('results/'+objective+'predicted_'+objective+'_inputs_guess='+str(num_guess)+'.csv', header=header)
    elif epoch > 0:
        input_parameter_trace.to_csv('results/'+objective+'predicted_'+objective+'_inputs_guess='+str(num_guess)+'.csv', mode='a', header=header)


def main():

    delete_previous_results(objective)
    model = torch.load("models/model.pt",map_location=device)

    # Load data
    features_scaling, e_scaling, V_scaling = getNormalization()
    features_train, R1_train, V_train, R2_train, e_train = getDataset(features_scaling, e_scaling, V_scaling, return_individual=True)
    
    e_train = e_scaling.unnormalize(e_train)
    e_train_matrix = assemble_e_matrix(e_train)
    train_objectives = loss_func(e_train_matrix, objective)
    _, best_idx = torch.topk(train_objectives, opt_num_guesses)

    print('Best performance in training dataset for objective '+objective+': '+str(torch.max(e_train[best_idx]).item()))

    features_best = features_train[best_idx,:]
    R1_best = R1_train[best_idx,:]
    R2_best = R2_train[best_idx,:]
    V_best = V_train[best_idx,:]
    R1_best_matrix = get_rotation_matrix(R1_best[:,0],R1_best[:,1],R1_best[:,2]).float()
    R1_best_6D = matrix_to_rotation_6d(R1_best_matrix).float()
    R2_best_matrix = get_rotation_matrix(R2_best[:,0],R2_best[:,1],R2_best[:,2]).float()
    R2_best_6D = matrix_to_rotation_6d(R2_best_matrix).float()
    input_best = torch.cat((features_best.to(device),R1_best_6D.to(device),V_best.to(device),R2_best_6D.to(device)),dim=1)


    for num_guess in range(opt_num_guesses):

        print('\nEpoch xxx | '+objective)
        predicted_parameters = input_best[num_guess,:].to(device)
        predicted_parameters.requires_grad = True
        optimizer = optim.Adam([predicted_parameters],lr=opt_lr)

        epoch = 0
        while epoch < opt_max_epochs:
            loss = get_loss(model, e_scaling, predicted_parameters, epoch)

            save_trace(predicted_parameters, features_scaling, V_scaling,epoch, num_guess)
            
            loss.backward()
            optimizer.step(lambda: get_loss(model, e_scaling, predicted_parameters, epoch))
            optimizer.zero_grad()
            epoch = epoch + 1
            if epoch % (opt_max_epochs/10) == 0:
                print("Epoch "+str(epoch)+": "+str(-loss.item()))

if __name__ == "__main__":
    main()