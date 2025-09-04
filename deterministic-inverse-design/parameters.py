import torch

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

################ Laptop########################
dataPath_input = 'data/dataset_input.csv'
dataPath_output = 'data/dataset_output.csv'

# define column names from data
features_names = ['relative_density','U1','U2','U3','lattice_type1','lattice_type2','lattice_type3','lattice_rep1','lattice_rep2','lattice_rep3']
R1_names = ['R1_theta','R1_rot_ax1','R1_rot_ax2']
V_names = ['V1','V2','V3']
R2_names = ['R2_theta','R2_rot_ax1','R2_rot_ax2']
all_names = features_names + R1_names + R2_names + V_names
C_ort_names = ['C11','C12','C13','C22','C23','C33','C44','C55','C66']
e_ort_names = ['e15','e24','e31','e32','e33']
C_names = ['C11','C12','C13','C14','C15','C16','C22','C23','C24','C25','C26','C33','C34','C35','C36','C44','C45','C46','C55','C56','C66']
e_names = ['e11','e12','e13','e14','e15','e16','e21','e22','e23','e24','e25','e26','e31','e32','e33','e34','e35','e36']
mu_names = ['mu11','mu22','mu33']
e_tilde_names = e_names

# define data type for correct scaling
features_types = ['continuous']*4 + ['categorical']*6
V_types = ['continuous']*3
# C_ort_types = ['continuous']*9
C_types = ['continuous']*21
e_types = ['continuous']*18
e_tilde_types = ['continuous']*18
C_ort_types = ['continuous']*9
e_ort_types = ['continuous']*5
mu_types = ['continuous']*3
d_types = ['continuous']*18

# define scaling strategies
features_scaling_strategy = 'min-max-1'
V_scaling_strategy = 'min-max-1'
C_ort_scaling_strategy = 'min-max-2'
e_ort_scaling_strategy = 'min-max-2'
e_tilde_scaling_strategy = 'min-max-2'
C_scaling_strategy = 'min-max-2'
C_hat_scaling_strategy = 'min-max-2'
e_scaling_strategy = 'min-max-2'
e_hat_scaling_strategy = 'min-max-2'
d_scaling_strategy = 'min-max-2'
mu_scaling_strategy = 'min-max-2'

# define train/test split
traintest_split = 0.99

# define batch size and numWorkers
batchSize = 8096
numWorkers = 0


train_epochs = 300
# F2_train_epochs = 200
model_arch = [1024,'leaky',1024,'leaky',1024,'leaky',1024,'leaky',1024,'leaky',1024,'leaky',1024]
learning_rate = 1e-3


# define loss function
lossFn = torch.nn.MSELoss()
shrinkage_scaler = 0.3
magnitude_scaler = 0.3
limit_stretch = True

opt_num_guesses = 150
opt_max_epochs = 10000
opt_lr = 0.01
use_dataset_guesses = True
