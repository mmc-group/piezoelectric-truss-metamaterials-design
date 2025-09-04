from math import *
import numpy as np
import torch
import torch.nn.functional as F
from parameters import *

Voigt_notation = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]

def rotate_elastic_constants(C, R):
    return tensor_to_Voigt(torch.einsum('...ia,...jb,...kc,...ld,...abcd->...ijkl',
                                               R,R,R,R,
                                               Voigt_to_tensor(C)))

def tensor_to_Voigt(C):
    batch_size = len(C)
    Voigt = torch.zeros((batch_size,6,6), dtype=float, device=device)
    for i in range(6):
        for j in range(6):
            k, l = Voigt_notation[i]
            m, n = Voigt_notation[j]
            Voigt[:,i,j] = C[:,k,l,m,n]
    return Voigt

def Voigt_to_tensor(C):
    batch_size = len(C)
    C_out = torch.zeros((batch_size,3,3,3,3), dtype=float, device=device)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    Voigt_i = full_index_to_Voigt_index(i,j)
                    Voigt_j = full_index_to_Voigt_index(k,l)
                    C_out[:,i,j,k,l] = C[:,Voigt_i,Voigt_j]
    return C_out

def full_index_to_Voigt_index(i, j):
    if i == j:
        return i
    return 6-i-j

def Voigt_to_Voigt_21(C):
    batch_size = len(C)
    C_out = torch.zeros((batch_size,21), dtype=float, device=device)
    iter = 0
    for i in range(6):
        for j in range (i,6):
            C_out[:,iter] = C[:,i,j]
            iter += 1
    return C_out

def e_Voigt_to_18(e):
    batch_size = len(e)
    e_out = torch.zeros((batch_size,18), dtype=float, device=device)
    iter = 0
    for i in range(3):
        for j in range (6):
            e_out[:,iter] = e[:,i,j]
            iter+=1
    return e_out

def get_rotation_matrix(theta, n1, n2):
    batch_size = len(theta)
    eps = 3e-5 # important to ensure n3-sqrt is not negative
    rot = torch.zeros((batch_size,3,3), dtype=float, device=device)
    # n3 = torch.sqrt(torch.nn.functional.relu(1. - torch.square(n1) - torch.square(n2)))
    n3 = torch.sqrt(1. - torch.square(n1) - torch.square(n2) + eps)
    # n_vec=torch.cat((n1.view(-1,1),n2.view(-1,1),n3.view(-1,1)),dim=1)
    # n_norm=torch.linalg.norm(n_vec,dim=1)
    # n1 = n1.clone()/torch.linalg.norm(n_vec)
    # n2 = n2.clone()/torch.linalg.norm(n_vec)
    # n3 = n3.clone()/torch.linalg.norm(n_vec)

    rot[:,0,0] = torch.square(n1)*(1-torch.cos(theta)) + torch.cos(theta)
    rot[:,0,1] = n1*n2*(1.-torch.cos(theta)) - n3*torch.sin(theta)
    rot[:,0,2] = n1*n3*(1.-torch.cos(theta)) + n2*torch.sin(theta)
    rot[:,1,0] = n2*n1*(1.-torch.cos(theta)) + n3*torch.sin(theta)
    rot[:,1,1] = torch.square(n2)*(1-torch.cos(theta)) + torch.cos(theta)
    rot[:,1,2] = n2*n3*(1.-torch.cos(theta)) - n1*torch.sin(theta)
    rot[:,2,0] = n3*n1*(1.-torch.cos(theta)) - n2*torch.sin(theta)
    rot[:,2,1] = n3*n2*(1.-torch.cos(theta)) + n1*torch.sin(theta)
    rot[:,2,2] = torch.square(n3)*(1.-torch.cos(theta)) + torch.cos(theta)
    # breakpoint()
    return rot

def get_rotation_matrix_np(theta, n1, n2):
    eps = 2e-6 # important to ensure n3-sqrt is not negative
    r = np.zeros((3,3))
    n3 = np.sqrt(1. - np.square(n1) - np.square(n2) + eps)
    r[0,0] = np.square(n1)*(1-np.cos(theta)) + np.cos(theta)
    r[0,1] = n1*n2*(1.-np.cos(theta)) - n3*np.sin(theta)
    r[0,2] = n1*n3*(1.-np.cos(theta)) + n2*np.sin(theta)
    r[1,0] = n2*n1*(1.-np.cos(theta)) + n3*np.sin(theta)
    r[1,1] = np.square(n2)*(1-np.cos(theta)) + np.cos(theta)
    r[1,2] = n2*n3*(1.-np.cos(theta)) - n1*np.sin(theta)
    r[2,0] = n3*n1*(1.-np.cos(theta)) - n2*np.sin(theta)
    r[2,1] = n3*n2*(1.-np.cos(theta)) + n1*np.sin(theta)
    r[2,2] = np.square(n3)*(1.-np.cos(theta)) + np.cos(theta)
    return r

def rotation_6d_to_matrix(d6: torch.Tensor) -> torch.Tensor:
    """
    Converts 6D rotation representation by Zhou et al. [1] to rotation matrix
    using Gram--Schmidt orthogonalisation per Section B of [1].
    Args:
        d6: 6D rotation representation, of size (*, 6)

    Returns:
        batch of rotation matrices of size (*, 3, 3)

    [1] Zhou, Y., Barnes, C., Lu, J., Yang, J., & Li, H.
    On the Continuity of Rotation Representations in Neural Networks.
    IEEE Conference on Computer Vision and Pattern Recognition, 2019.
    Retrieved from http://arxiv.org/abs/1812.07035
    """

    a1, a2 = d6[..., :3], d6[..., 3:]
    b1 = F.normalize(a1, dim=-1)
    b2 = a2 - (b1 * a2).sum(-1, keepdim=True) * b1
    b2 = F.normalize(b2, dim=-1)
    b3 = torch.cross(b1, b2, dim=-1)
    return torch.stack((b1, b2, b3), dim=-2)

def rot6DToAngleAxis(rot_6d):
    rot_matrix = rotation_6d_to_matrix(rot_6d)
    x = rot_matrix[:,1,2]-rot_matrix[:,2,1]
    y = rot_matrix[:,2,0]-rot_matrix[:,0,2]
    z = rot_matrix[:,0,1]-rot_matrix[:,1,0]
    r = torch.sqrt(torch.pow(x,2.)+torch.pow(y,2.)+torch.pow(z,2.))
    t = rot_matrix[:,0,0]+rot_matrix[:,1,1]+rot_matrix[:,2,2]
    theta = torch.atan2(r,t-1.)

    # We below enforce a positive e3-component of the angle-axis representation by mirroring the rotation axis w.r.t. the origin in the case e3 is indeed negative.
    # This is necessary as we construct e3 by e3=\sqrt(1.-e1^2-e2^2) (>= 0 for any e1 and e2).
    for i in range(len(z)):
        if (z[i] < 0):
            theta[i] = -theta[i]
            x[i] = -x[i]
            y[i] = -y[i]
    return torch.stack((-theta,torch.div(x,r),torch.div(y,r)),dim=1)

def matrix_to_rotation_6d(rot):
    return torch.flatten(torch.split(rot,[2,1],dim=1)[0],start_dim=1, end_dim=2).float()

# direct piezo tensor rotation
def direct_rotate_e(e,rot,orthotropic=False,method=None):

    R = get_rotation_matrix(rot[:,0],rot[:,1],rot[:,2]).float()
    batch_size = len(e)
    e_Voigt = torch.zeros((batch_size,3,6), device=device)

    if orthotropic:
        e_Voigt[:,0,4] = e[:,0]
        e_Voigt[:,1,3] = e[:,1]
        e_Voigt[:,2,0] = e[:,2]
        e_Voigt[:,2,1] = e[:,3]
        e_Voigt[:,2,2] = e[:,4]

    else:
        e_Voigt[:,0,0] = e[:,0]
        e_Voigt[:,0,1] = e[:,1]
        e_Voigt[:,0,2] = e[:,2]
        e_Voigt[:,0,3] = e[:,3]
        e_Voigt[:,0,4] = e[:,4]
        e_Voigt[:,0,5] = e[:,5]

        e_Voigt[:,1,0] = e[:,6]
        e_Voigt[:,1,1] = e[:,7]
        e_Voigt[:,1,2] = e[:,8]
        e_Voigt[:,1,3] = e[:,9]
        e_Voigt[:,1,4] = e[:,10]
        e_Voigt[:,1,5] = e[:,11]

        e_Voigt[:,2,0] = e[:,12]
        e_Voigt[:,2,1] = e[:,13]
        e_Voigt[:,2,2] = e[:,14]
        e_Voigt[:,2,3] = e[:,15]
        e_Voigt[:,2,4] = e[:,16]
        e_Voigt[:,2,5] = e[:,17]



    # breakpoint()
    K0 = torch.mul(R,R)
    K1 = torch.zeros((batch_size,3,3), device=device)
    K2 = torch.zeros((batch_size,3,3), device=device)
    K3 = torch.zeros((batch_size,3,3), device=device)
    K = torch.zeros((batch_size,6,6), device=device)

    K1[:,0,0] = R[:,0,1]*R[:,0,2]
    K1[:,0,1] = R[:,0,2]*R[:,0,0]
    K1[:,0,2] = R[:,0,0]*R[:,0,1]
    K1[:,1,0] = R[:,1,1]*R[:,1,2]
    K1[:,1,1] = R[:,1,2]*R[:,1,0]
    K1[:,1,2] = R[:,1,0]*R[:,1,1]
    K1[:,2,0] = R[:,2,1]*R[:,2,2]
    K1[:,2,1] = R[:,2,2]*R[:,2,0]
    K1[:,2,2] = R[:,2,0]*R[:,2,1]

    K2[:,0,0] = R[:,1,0]*R[:,2,0]
    K2[:,0,1] = R[:,1,1]*R[:,2,1]
    K2[:,0,2] = R[:,1,2]*R[:,2,2]
    K2[:,1,0] = R[:,2,0]*R[:,0,0]
    K2[:,1,1] = R[:,2,1]*R[:,0,1]
    K2[:,1,2] = R[:,2,2]*R[:,0,2]
    K2[:,2,0] = R[:,0,0]*R[:,1,0]
    K2[:,2,1] = R[:,0,1]*R[:,1,1]
    K2[:,2,2] = R[:,0,2]*R[:,1,2]

    K3[:,0,0] = R[:,1,1]*R[:,2,2]+R[:,1,2]*R[:,2,1]
    K3[:,0,1] = R[:,1,2]*R[:,2,0]+R[:,1,0]*R[:,2,2]
    K3[:,0,2] = R[:,1,0]*R[:,2,1]+R[:,1,1]*R[:,2,0]
    K3[:,1,0] = R[:,2,1]*R[:,0,2]+R[:,2,2]*R[:,0,1]
    K3[:,1,1] = R[:,2,2]*R[:,0,0]+R[:,2,0]*R[:,0,2]
    K3[:,1,2] = R[:,2,0]*R[:,0,1]+R[:,2,1]*R[:,0,0]
    K3[:,2,0] = R[:,0,1]*R[:,1,2]+R[:,0,2]*R[:,1,1]
    K3[:,2,1] = R[:,0,2]*R[:,1,0]+R[:,0,0]*R[:,1,2]
    K3[:,2,2] = R[:,0,0]*R[:,1,1]+R[:,0,1]*R[:,1,0]

    K1 = 2*K1
    K_top = torch.cat((K0,K1),2)
    K_bot = torch.cat((K2,K3),2)
    K = torch.cat((K_top,K_bot),1).float()
    # breakpoint()
    breakpoint()
    return e_Voigt_to_18(torch.matmul(torch.matmul(R,e_Voigt),torch.transpose(K,1,2))).float()
    # return e_Voigt_to_18(torch.matmul(torch.matmul(R,e_Voigt),K)).float()

# direct stiffness rotation in Voigt notation (for details see Bower's 'Applied Mechanics of Solids' [Chapter 3])
def direct_rotate(C,rot,orthotropic=False,method=None):

    if method == '6D':
        R = rotation_6d_to_matrix(rot)
    else:
        R = get_rotation_matrix(rot[:,0],rot[:,1],rot[:,2])
    batch_size = len(C)
    Voigt = torch.zeros((batch_size,6,6), device=device)

    if orthotropic:
        Voigt[:,0,0] = C[:,0]
        Voigt[:,0,1] = C[:,1]
        Voigt[:,1,0] = C[:,1]
        Voigt[:,0,2] = C[:,2]
        Voigt[:,2,0] = C[:,2]
        Voigt[:,1,1] = C[:,3]
        Voigt[:,1,2] = C[:,4]
        Voigt[:,2,1] = C[:,4]
        Voigt[:,2,2] = C[:,5]
        Voigt[:,3,3] = C[:,6]
        Voigt[:,4,4] = C[:,7]
        Voigt[:,5,5] = C[:,8]
    else:
        Voigt[:,0,0] = C[:,0]
        Voigt[:,0,1] = C[:,1]
        Voigt[:,1,0] = C[:,1]
        Voigt[:,0,2] = C[:,2]
        Voigt[:,2,0] = C[:,2]
        Voigt[:,0,3] = C[:,3]
        Voigt[:,3,0] = C[:,3]
        Voigt[:,0,4] = C[:,4]
        Voigt[:,4,0] = C[:,4]
        Voigt[:,0,5] = C[:,5]
        Voigt[:,5,0] = C[:,5]

        Voigt[:,1,1] = C[:,6]
        Voigt[:,1,2] = C[:,7]
        Voigt[:,2,1] = C[:,7]
        Voigt[:,1,3] = C[:,8]
        Voigt[:,3,1] = C[:,8]
        Voigt[:,1,4] = C[:,9]
        Voigt[:,4,1] = C[:,9]
        Voigt[:,1,5] = C[:,10]
        Voigt[:,5,1] = C[:,10]

        Voigt[:,2,2] = C[:,11]
        Voigt[:,2,3] = C[:,12]
        Voigt[:,3,2] = C[:,12]
        Voigt[:,2,4] = C[:,13]
        Voigt[:,4,2] = C[:,13]
        Voigt[:,2,5] = C[:,14]
        Voigt[:,5,2] = C[:,14]

        Voigt[:,3,3] = C[:,15]
        Voigt[:,3,4] = C[:,16]
        Voigt[:,4,3] = C[:,16]
        Voigt[:,3,5] = C[:,17]
        Voigt[:,5,3] = C[:,17]

        Voigt[:,4,4] = C[:,18]
        Voigt[:,4,5] = C[:,19]
        Voigt[:,5,4] = C[:,19]

        Voigt[:,5,5] = C[:,20]

    K0 = torch.mul(R,R)
    K1 = torch.zeros((batch_size,3,3), device=device)
    K2 = torch.zeros((batch_size,3,3), device=device)
    K3 = torch.zeros((batch_size,3,3), device=device)
    K = torch.zeros((batch_size,6,6), device=device)

    K1[:,0,0] = R[:,0,1]*R[:,0,2]
    K1[:,0,1] = R[:,0,2]*R[:,0,0]
    K1[:,0,2] = R[:,0,0]*R[:,0,1]
    K1[:,1,0] = R[:,1,1]*R[:,1,2]
    K1[:,1,1] = R[:,1,2]*R[:,1,0]
    K1[:,1,2] = R[:,1,0]*R[:,1,1]
    K1[:,2,0] = R[:,2,1]*R[:,2,2]
    K1[:,2,1] = R[:,2,2]*R[:,2,0]
    K1[:,2,2] = R[:,2,0]*R[:,2,1]

    K2[:,0,0] = R[:,1,0]*R[:,2,0]
    K2[:,0,1] = R[:,1,1]*R[:,2,1]
    K2[:,0,2] = R[:,1,2]*R[:,2,2]
    K2[:,1,0] = R[:,2,0]*R[:,0,0]
    K2[:,1,1] = R[:,2,1]*R[:,0,1]
    K2[:,1,2] = R[:,2,2]*R[:,0,2]
    K2[:,2,0] = R[:,0,0]*R[:,1,0]
    K2[:,2,1] = R[:,0,1]*R[:,1,1]
    K2[:,2,2] = R[:,0,2]*R[:,1,2]

    K3[:,0,0] = R[:,1,1]*R[:,2,2]+R[:,1,2]*R[:,2,1]
    K3[:,0,1] = R[:,1,2]*R[:,2,0]+R[:,1,0]*R[:,2,2]
    K3[:,0,2] = R[:,1,0]*R[:,2,1]+R[:,1,1]*R[:,2,0]
    K3[:,1,0] = R[:,2,1]*R[:,0,2]+R[:,2,2]*R[:,0,1]
    K3[:,1,1] = R[:,2,2]*R[:,0,0]+R[:,2,0]*R[:,0,2]
    K3[:,1,2] = R[:,2,0]*R[:,0,1]+R[:,2,1]*R[:,0,0]
    K3[:,2,0] = R[:,0,1]*R[:,1,2]+R[:,0,2]*R[:,1,1]
    K3[:,2,1] = R[:,0,2]*R[:,1,0]+R[:,0,0]*R[:,1,2]
    K3[:,2,2] = R[:,0,0]*R[:,1,1]+R[:,0,1]*R[:,1,0]

    K1 = 2*K1
    K_top = torch.cat((K0,K1),2)
    K_bot = torch.cat((K2,K3),2)
    K = torch.cat((K_top,K_bot),1).float()

    return Voigt_to_Voigt_21(torch.matmul(torch.matmul(K,Voigt),torch.transpose(K,1,2))).float()
