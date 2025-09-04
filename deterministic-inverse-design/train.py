import pathlib
import torch
from torch.utils.data import DataLoader
from parameters import *
from src.loadDataset import *
from src.model_utils import *
from src.errorAnalysis import computeR2

if __name__ == '__main__':
    torch.cuda.empty_cache()
    np.random.seed(1234)
    torch.manual_seed(1234)
    loss_best = torch.inf
    # create directories
    pathlib.Path('models').mkdir(exist_ok=True)
    pathlib.Path('training').mkdir(exist_ok=True)
    pathlib.Path('training/history').mkdir(exist_ok=True)

    # load and preprocess data
    features_scaling, e_scaling, V_scaling = getNormalization()
    train_set, test_set = getDataset(features_scaling, e_scaling, V_scaling)
    train_data_loader = DataLoader(dataset=train_set, num_workers=numWorkers, batch_size=batchSize)
    test_data_loader = DataLoader(dataset=test_set, num_workers=numWorkers, batch_size=len(test_set))
    # Note: for test, batch_size=len(test_set) so that we load the entire test set at once
    features_test, R1_test, V_test, R2_test, e_test = next(iter(test_data_loader))
    print('\n-------------------------------------')

    # initialize first forward model
    model = createNN(46,model_arch,18).to(device)
    # set up optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    train_history, test_history, R2_history = [],[],[]
    # training
    for epoch_iter in range(train_epochs):
        train_loss = 0.
        for iteration, batch in enumerate(train_data_loader,0):
            # get batch
            features_train, R1_train, V_train, R2_train, e_train = batch[0].to(device), batch[1].to(device),\
                            batch[2].to(device), batch[3].to(device), batch[4].to(device)

            R1_train_matrix = get_rotation_matrix(R1_train[:,0],R1_train[:,1],R1_train[:,2]).float()
            R1_train_6D = matrix_to_rotation_6d(R1_train_matrix).float()
            R2_train_matrix = get_rotation_matrix(R2_train[:,0],R2_train[:,1],R2_train[:,2]).float()
            R2_train_6D = matrix_to_rotation_6d(R2_train_matrix).float()
            # set train mode
            model.train()
            # forward pass F1
            e_train_pred = model(torch.cat((features_train,R1_train_6D,V_train,R2_train_6D),dim=1))
            loss = lossFn(e_train_pred,e_train)
            # optimize
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            # store (batch) training loss
            train_loss = loss.item()

        features_test, V_test, R1_test, R2_test, e_test = features_test.to(device), V_test.to(device), R1_test.to(device), R2_test.to(device), e_test.to(device)
        R1_test_matrix = get_rotation_matrix(R1_test[:,0],R1_test[:,1],R1_test[:,2]).float()
        R1_test_6D = matrix_to_rotation_6d(R1_test_matrix).float()
        R2_test_matrix = get_rotation_matrix(R2_test[:,0],R2_test[:,1],R2_test[:,2]).float()
        R2_test_6D = matrix_to_rotation_6d(R2_test_matrix).float()
        e_test_pred = model(torch.cat((features_test, R1_test_6D, V_test, R2_test_6D), dim=1))

        test_loss = lossFn(e_test_pred,e_test).item()

        print("| {}:{}/{} | F1_EpochTrainLoss: {:.2e} | F1_EpochTestLoss: {:.2e} | F1-R2 e: {:.2f}".format("F1", epoch_iter, train_epochs,
                    train_loss, test_loss,
                    torch.mean(computeR2(e_test_pred, e_test)).item()))
        train_history.append(train_loss)
        test_history.append(test_loss)
        R2_history.append(torch.mean(computeR2(e_test_pred, e_test)).item())
        # save model
        if test_loss < loss_best:
            loss_best = test_loss
            torch.save(model,"models/model.pt")
        # export loss history
        exportList('training/history/train_loss_history',train_history)
        exportList('training/history/test_loss_history',test_history)
        exportList('training/history/R2_history',R2_history)

    print('\n-------------------------------------')

    torch.cuda.empty_cache()
    print('Finished.')

    