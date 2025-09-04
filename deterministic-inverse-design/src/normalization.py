import torch
import torch.nn.functional as F
import numpy as np
from parameters import *

class Normalization:
    def __init__(self,data,dataType,strategy):
        self.mu = torch.mean(data,dim=0)
        self.std = torch.std(data,dim=0)
        self.min = torch.min(data,dim=0)[0]
        self.max = torch.max(data,dim=0)[0]
        self.globalmin = torch.min(data)
        self.globalmax = torch.max(data)
        self.dataType = dataType
        self.cols = data.size()[1]
        self.strategy = strategy

    def normalize(self, data):
        list_index_cat = []
        temp = torch.zeros(data.shape,device=data.device)
        for i in range(0, self.cols):
            if self.dataType[i] == 'continuous':

                if(self.strategy == 'min-max-1'):
                    #scale to [0,1]
                    temp[:,i] = torch.div(data[:,i]-self.min[i], self.max[i]-self.min[i])

                elif(self.strategy == 'global-min-max-1'):
                    #scale to [-1,1] based on min max of full dataset
                    temp[:,i] = torch.div(data[:,i]-self.globalmin, self.globalmax-self.globalmin)

                elif(self.strategy == 'min-max-2'):
                    #scale to [-1,1]
                    temp[:,i] = 2.*torch.div(data[:,i]-self.min[i], self.max[i]-self.min[i])-1.

                elif(self.strategy == 'global-min-max-2'):
                    #scale to [-1,1] based on min max of full dataset
                    temp[:,i] = 2.*torch.div(data[:,i]-self.globalmin, self.globalmax-self.globalmin)-1.

                elif(self.strategy == 'mean-std'):
                    #scale s.t. mean=0, std=1
                    temp[:,i] = torch.div(data[:,i]-self.mu[i], self.std[i])

                elif (self.strategy == 'none'):
                    temp[:,i] = data[:,i]

                else:
                    raise ValueError('Incorrect normalization strategy')

            elif self.dataType[i] == 'categorical':
                num_classes = self.max[i]
                #convert categorical features into binaries and append at the end of feature tensor
                temp = torch.cat((temp,F.one_hot(data[:,i].to(torch.int64),num_classes=int(num_classes.item()+1))),dim=1)
                list_index_cat = np.append(list_index_cat,i)

            else:
                raise ValueError("Data type must be either continuous or categorical")

        # delete original (not one-hot encoded) categorical features
        j = 0
        for i in np.array(list_index_cat, dtype=np.int64):
            temp = torch.cat([temp[:,0:i+j], temp[:,i+1+j:]],dim=1)
            j -= 1

        return temp

    def unnormalize(self, data):
        temp = torch.zeros(data.shape,device=data.device)
        for i in range(0, self.cols):
            if self.dataType[i] == 'continuous':

                if(self.strategy == 'min-max-1'):
                    temp[:,i] = torch.mul(data[:,i], self.max[i]-self.min[i]) +self.min[i]

                elif(self.strategy == 'global-min-max-1'):
                    temp[:,i] = torch.mul(data[:,i], self.globalmax-self.globalmin) +self.globalmin

                elif(self.strategy == 'min-max-2'):
                    temp[:,i] = torch.mul(0.5*data[:,i]+0.5, self.max[i]-self.min[i]) +self.min[i]

                elif(self.strategy == 'global-min-max-2'):
                    temp[:,i] = torch.mul(0.5*data[:,i]+0.5, self.globalmax-self.globalmin) +self.globalmin

                elif(self.strategy == 'mean-std'):
                    temp[:,i] = torch.mul(data[:,i], self.std[i]) + self.mu[i]

                elif (self.strategy == 'none'):
                    temp[:,i] = data[:,i]

                else:
                    raise ValueError('Incorrect normalization strategy')

            elif self.dataType[i] == 'categorical':
                temp[:,i] = data[:,i]

            else:
                raise ValueError("Data type must be either continuous or categorical")
        return temp

# convert one-hot representation back to categorical integers
def decodeOneHot(data):
    type1,type2,type3,rep1,rep2,rep3 = torch.split(data,[7,7,7,2,2,2],dim=1)
    # we increment the repetition to convert from binary to [1,2] as defined in the publication
    type1,type2,type3,rep1,rep2,rep3 = torch.argmax(type1, dim=1),torch.argmax(type2,dim=1),torch.argmax(type3,dim=1),torch.argmax(rep1,dim=1)+1,torch.argmax(rep2,dim=1)+1,torch.argmax(rep3, dim=1)+1

    types = torch.stack((type1,type2,type3),dim=1)
    reps = torch.stack((rep1,rep2,rep3),dim=1)

    # sort by lattice number
    sorted_types, indices = torch.sort(types)
    sorted_reps = smart_sort(reps, indices)

    # sort by repetitions if lattice numbers are equal
    for i in range(data.size()[0]):
        if sorted_types[i,0] == sorted_types[i,1] and sorted_types[i,1] == sorted_types[i,2]:
            sorted_reps[i,:] = torch.sort(sorted_reps[i,:])[0]
        elif sorted_types[i,0] == sorted_types[i,1]:
            sorted_reps[i,0:2] = torch.sort(sorted_reps[i,0:2])[0]
        elif sorted_types[i,1] == sorted_types[i,2]:
            sorted_reps[i,1:3] = torch.sort(sorted_reps[i,1:3])[0]
        else:
            pass

    return torch.cat((sorted_types,sorted_reps), dim=1)

def smart_sort(x, permutation):
    d1, d2 = x.size()
    ret = x[
        torch.arange(d1).unsqueeze(1).repeat((1, d2)).flatten(),
        permutation.flatten()
    ].view(d1, d2)
    return ret
