clc; clear all; close all
addpath("..\data")
warning ('off','all')
columns_prop = {'C11','C12','C13','C14','C15','C16','C21','C22','C23', ...
    'C24','C25','C26','C31', 'C32', 'C33','C34','C35','C36',...
    'C41', 'C42', 'C43', 'C44','C45','C46', ...
    'C51', 'C52', 'C53', 'C54', 'C55','C56',...
    'C61', 'C62', 'C63', 'C64', 'C65', 'C66',...
    'e11','e12','e13','e14','e15','e16','e21', ...
    'e22','e23','e24','e25','e26','e31','e32','e33','e34','e35', ...
    'e36', ...
    'mu11','mu12','mu13',...
    'mu21','mu22','mu23',...
    'mu31','mu32','mu33'};
prop_table = array2table(zeros(0, numel(columns_prop)), 'VariableNames', columns_prop);
%%% DSI or DSII for data generation for dataset I or II
dataset="DSI";
if dataset=="DSI"
    load input_DS1.mat
    range=1:size(refdata,1);
    train_dens=refdata(range,1);
    train_Rot1=refdata(range,2:4);train_Rot2=refdata(range,5:7);
    train_stretch1=refdata(range,8:10);train_stretch2=refdata(range,11:13);
    train_lat=refdata(range,14:2:19);train_tess=refdata(range,15:2:19);
elseif dataset=="DSII"
    load nodes_DS2.mat
    load connect_DS2.mat
    num_trusses=length(nodes_full);
    range=1:num_trusses;
    rho=0.1;  %Density of the lattice
end
P_mat=zeros(length(range),63);

%%% material properties
E=60.606e9; % Young's Modulus
nu=0.3;  %Poisson's Ratio
mu33=1433.6*8.854e-12;    %permittivity
% piezoelectric coefficients' matrix (stress-charge form)
e_base=[0	0	0	0	17.0345	0
    0	0	0	17.0345	0	0
    -6.62281	-6.62281	23.2403	0	0	0];
for tr=1:length(range)
    if dataset == "DSI"
        rho=train_dens(tr,:);  %Density of the lattice
        THETA.lat=train_lat(tr,:);
        THETA.tess=train_tess(tr,:);
        THETA.stretch1=train_stretch1(tr,:);THETA.stretch2=train_stretch2(tr,:);
        THETA.Rot1=train_Rot1(tr,:);THETA.Rot2=train_Rot2(tr,:);
        [L1,L2,L3]=lattice_vectors();
        [Coords,Connectivity,Lattice_vectors]=Generate_Lattice(THETA,L1,L2,L3);
        L1=Lattice_vectors(1,:);
        L2=Lattice_vectors(2,:);
        L3=Lattice_vectors(3,:);
    else
        Coords = nodes_full{tr};
        Connectivity = connect_full{tr};
        [Coords]=ScaleToUnit(Coords,Connectivity);
        [L1,L2,L3]=lattice_vectors();
    end
    %%% homogenization for effective properties
    [Relative_density,C,e,mu,d]=homogenization(rho,Coords,Connectivity,L1,L2,L3,E,nu,mu33,e_base);
    P_mat(tr,1:36)=reshape(C',[36,1]);
    P_mat(tr,37:54)=reshape(e',[18,1]);
    P_mat(tr,55:63)=reshape(mu',[9,1]);
    %%% conversion to compliance (strain-charge) form
    % S=inv(C);
    % d=e*S;

    %%% hydrostatic piezoelectric coefficient
    % fomhs(tr)=e(3,1)+e(3,2)+e(3,3);

    %%% plotting effective elasticity and piezoelectric surfaces (disabled during dataset generation)
    % effective_piezo_surface(e)
    % effective_surface_plot(C)

    %%% visualize the unit cell
    % truss_plotting(Coords, Connectivity)
    %%% storing properties
    new_row = array2table(P_mat(tr,:), 'VariableNames', columns_prop);
    prop_table = [prop_table; new_row];

    disp(tr)
end
%%% saving the dataset in a file
% writetable(prop_table,'dataset.csv')
