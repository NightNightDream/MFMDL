import torch
import pandas as pd
import numpy as np
import random

from torch.nn import Linear
from torch_geometric.nn import GCNConv
from torch.utils.data import DataLoader 
from torch.utils.data import Dataset
from scipy.sparse import coo_matrix
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_curve, auc

class MyGCN(torch.nn.Module):
    def __init__(self, adj, gene_num):
        super().__init__()
        self.conv1 = MyGCNConv(1, 1, adj, True)
        self.linear1 = Linear(gene_num, 256)
        self.linear2 = Linear(256, 128)
        self.linear3 = Linear(128, 32)

        self.linear_TILs = Linear(28, 16)

        self.linear_PW = Linear(50, 32)

        self.linear_pre = Linear(80, 1)

        self.dropout = torch.nn.Dropout(p=drop_rate)


    
    def forward(self, nodes_feature, TILs_feature, pws_feature, predict=False):
        h = nodes_feature
        h = torch.nn.functional.normalize(h, dim=1)
        # print(h.shape)
        # print(torch.squeeze(h))
        h = self.conv1(h)
        batch = torch.squeeze(h)
        # print(batch)
        out1 = self.dropout(self.linear1(batch))
        out1 = self.dropout(self.linear2(out1))
        out1 = self.dropout(self.linear3(out1))
        out1 = torch.relu(out1)
        # print(out1)
        if predict:
            out1 = torch.stack([out1])

        out2 = self.linear_TILs(TILs_feature)
        out2 = torch.relu(out2)
        #print(out2)
        out3 = self.linear_PW(pws_feature)
        out3 = torch.relu(out3)
        
        out = self.linear_pre(torch.cat((out1, out2, out3), dim=1))
        out = torch.sigmoid(out)

        # print(torch.cat((out1, out2, out3), dim=1))
        # print(torch.squeeze(out))
        return torch.squeeze(out)


class MyGCNConv(torch.nn.Module):
    def __init__(self, in_features, out_features, adj, bias=True):
        super().__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.adj = adj
        if bias:
            self.bias = torch.nn.Parameter(torch.FloatTensor(out_features))
        else:
            self.register_parameter('bias', None)

    def forward(self, node_features):
        batch = []
        for i in range(node_features.shape[0]):
            out = torch.matmul(self.adj, node_features[i])
            if self.bias is not None:
                batch.append(out + self.bias)
            else:
                batch.append(out)
        return torch.stack(batch)

class MyDataset(Dataset):
    def __init__(self, nodes_feature, TILs_feature, pws_feature, labels):
        self.nodes_feature = nodes_feature
        self.TILs_feature = TILs_feature
        self.pws_feature = pws_feature
        self.labels = labels
    
    def __getitem__(self, index):
        return self.nodes_feature[index], self.labels[index], self.TILs_feature[index], self.pws_feature[index]

    def __len__(self):
        return len(self.nodes_feature)



regulate_network = pd.read_csv("./data/Hand_Data/regulate_network.txt", sep=",")

def get_interacter_genes(dataset_name):
    # network_genes = list(set(regulate_network['Gene1'].tolist()) | set(regulate_network["Gene2"].tolist()))
    network_genes = list(pd.read_csv("./data/Hand_Data/selected_genes.csv", sep="\t")['Gene_Symbol'].tolist())
    expr_genes = pd.read_csv("./data/Sour_Data/" + dataset_name + "/" + dataset_name + "_norm3.txt", sep="\t")
    return list(set(network_genes) & set(expr_genes.loc[:,'gene_id'].tolist()))



def get_adj_matrix():
    gene_dict = {}
    for i in range(len(selected_genes)):
        gene_dict[selected_genes[i]] = i
    
    edge_index = [[],[]]
    edge_weight = []
    for i in range(len(regulate_network)):
        gene1 = regulate_network.iloc[i, 0]
        gene2 = regulate_network.iloc[i, 1]
        corVal = regulate_network.iloc[i, 2]
        if gene_dict.__contains__(gene1) and gene_dict.__contains__(gene2):
            edge_index[0].append(gene_dict[gene1])
            edge_index[1].append(gene_dict[gene2])
            edge_weight.append(corVal)
    edge_index = torch.tensor(edge_index)
    edge_weight = torch.tensor(edge_weight, dtype=torch.float32)

    dim = len(gene_dict)
    coo_rows = np.array(edge_index[0])
    coo_cols = np.array(edge_index[1])
    coo_data = np.array(edge_weight)
    coo_edge = np.ones(len(edge_weight))

    adj = coo_matrix((coo_edge, (coo_rows, coo_cols)), shape=(dim, dim)).todense() + np.eye(dim)
    D_sqrt = np.diag(np.power(np.array(adj.sum(1)), -0.5).flatten())
    temp = np.matmul(D_sqrt, coo_matrix((coo_data, (coo_rows, coo_cols)), shape=(dim, dim)).todense() + np.eye(dim))
    # print(np.matmul(temp, D_sqrt)[3])
    return np.matmul(temp, D_sqrt)
    

class DataHandler:
    dataset_response_labels = []
    dataset_nodes_feature = []
    dataset_TILs_feature = []
    dataset_pws_feature = []

    @staticmethod
    def get_dataloader(dataset_numbers, train_idx, test_idx):
        
        if len(DataHandler.dataset_response_labels) == 0:
            # select top 50
            pw_df = pd.read_csv("./data/Result/pw_feature.txt", sep="\t")
            pws_name = pw_df.iloc[0:50,0].tolist()

            for i in range(len(dataset_numbers)):
                path = "./data/Sour_Data/" + dataset_numbers[i] + "/" + dataset_numbers[i]
                response_df = pd.read_csv(path + "_patient.txt", sep="\t")
                sample_names = response_df['Patient'].tolist() 
                if dataset_numbers[i] == "GSE115821":
                    for j in range(len(sample_names)):
                        sample_names[j] = sample_names[j].replace("-",".")
                DataHandler.dataset_response_labels.extend(response_df['flag'].tolist())

                gene_expr = pd.read_csv(path + "_norm3.txt", sep="\t")
                gene_expr = gene_expr.set_index('gene_id')
                gene_expr_temp = gene_expr.loc[selected_genes,]
                TILs_df = pd.read_csv(path + "_TILs.txt", sep="\t")
                for patient in sample_names:
                    temp = []
                    for expr in gene_expr_temp[patient].tolist():
                        temp.append(torch.tensor([expr], dtype=torch.float32))
                    
                    DataHandler.dataset_nodes_feature.append(torch.stack(temp))
                    DataHandler.dataset_TILs_feature.append(torch.tensor(TILs_df[patient], dtype=torch.float32))
                
                

                # pathway features
                ssgsea_df = pd.read_csv(path + "_ssgsea.txt", sep="\t")
                ssgsea_df = ssgsea_df.loc[pws_name,]
                for patinet in sample_names:
                    DataHandler.dataset_pws_feature.append(ssgsea_df.loc[:,patinet].tolist())


        dataset_response_labels = torch.tensor(DataHandler.dataset_response_labels, dtype=torch.float32)
        dataset_nodes_feature = torch.stack(DataHandler.dataset_nodes_feature)
        dataset_TILs_feature = torch.stack(DataHandler.dataset_TILs_feature)
        dataset_pws_feature = torch.tensor(DataHandler.dataset_pws_feature, dtype=torch.float32)

        split = int(len(train_idx) * 0.8)

        train_dataset = MyDataset(dataset_nodes_feature[train_idx[:split]],
                                dataset_TILs_feature[train_idx[:split]],
                                dataset_pws_feature[train_idx[:split]],
                                    dataset_response_labels[train_idx[:split]])
        train_dataloader = DataLoader(dataset=train_dataset, shuffle=False, batch_size=16, drop_last=False)

        valid_dataset = MyDataset(dataset_nodes_feature[train_idx[split:]],
                                dataset_TILs_feature[train_idx[split:]],
                                dataset_pws_feature[train_idx[split:]],
                                    dataset_response_labels[train_idx[split:]])
        valid_dataloader = DataLoader(dataset=valid_dataset, shuffle=False, batch_size=16, drop_last=False)

        test_dataset = MyDataset(dataset_nodes_feature[test_idx],
                                dataset_TILs_feature[test_idx],
                                dataset_pws_feature[test_idx],
                                    dataset_response_labels[test_idx])
        test_dataloader = DataLoader(dataset=test_dataset, shuffle=False, batch_size=16, drop_last=False)

        return train_dataloader, valid_dataloader, test_dataloader



def calculate_accuracy(model, dataloader, is_testdataset, print_flag = False):
    if is_testdataset:
        accuracy = 0
        for batch_idx, data in enumerate(dataloader, 0):
            nodes_feature, labels, TILs_feature, pws_feature = data
            output = model(nodes_feature, TILs_feature, pws_feature, True)
            if output > 0.5:
                predict = 1
                test_predicts.append(1)
            else:
                predict = 0
                test_predicts.append(0)
            test_proba.append(output.item())
            lable = 1
            if labels[0] == 1:
                test_labels.append(1)
            else:
                test_labels.append(0)
                lable = 0
            
        if print_flag:
            print("Pred: {} Label: {}".format(predict, lable))
        
        
    else:
        predicts = []
        accuracy = 0
        for batch_idx, data in enumerate(dataloader, 0):
            nodes_feature, labels, TILs_feature, pws_feature = data
            outputs = model(nodes_feature, TILs_feature, pws_feature, False)
            for i in range(len(outputs)):
                if outputs[i] > 0.5:
                    predict = 1
                else:
                    predict = 0
                predicts.append(predict)
                if predict == labels[i]:
                    accuracy += 1
        if print_flag:
            print(predicts)
    return accuracy / dataloader.dataset.__len__()


def get_valid_loss(model, cross_loss, valid_dl):
    val_loss = 0
    for batch_idx, data in enumerate(valid_dl, 0):
        nodes_feature, labels, TILs_feature, pws_feature = data
        model.eval()
        output = model(nodes_feature, TILs_feature, pws_feature)
        model.train()
        val_loss += cross_loss.item()
    # print(val_loss)
    return val_loss


# dataset_name = "Hugo"
# dataset_name = "Lee"
# dataset_name = "Liu"
# dataset_name = "Gide"
dataset_name = "Kim"

# hyperparameter
val_loss_threshold = 0.015
loss_threshold = 0.75
learning_rate = 0.0005
epoch_num = 200
drop_rate = 0.3


selected_genes = get_interacter_genes(dataset_name)
test_predicts = []
test_proba = []
test_labels = []
def run():
    torch.manual_seed(200)
    datasets = [dataset_name]
    path = "./data/Sour_Data/" + dataset_name + "/" + dataset_name
    loo = LeaveOneOut()
    for train_idx, test_idx in loo.split(pd.read_csv(path + "_patient.txt", sep="\t")['flag'].tolist()):
        flag = 1
        while flag == 1:
            for epoch in range(epoch_num):
                if epoch == 0:  # reset the model
                    random.shuffle(train_idx)
                    train_dl, valid_dl, test_dl = DataHandler.get_dataloader(datasets, train_idx, test_idx)
                    model = MyGCN(torch.tensor(get_adj_matrix(), dtype=torch.float), len(selected_genes))
                    criterion = torch.nn.BCELoss()
                    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, eps=1e-9)
                random.shuffle(train_idx)
                train_dl, valid_dl, test_dl = DataHandler.get_dataloader(datasets, train_idx, test_idx)
                train_loss = 0.0
                for batch_idx, data in enumerate(train_dl, 0):
                    nodes_feature, labels, TILs_feature, pws_feature = data
                    optimizer.zero_grad()
                    output = model(nodes_feature, TILs_feature, pws_feature)
                    try:
                        cross_loss = criterion(output, labels)
                    except:
                        print(output)
                    # l1_regularization, l2_regularization = torch.tensor(0, dtype=torch.float32), torch.tensor(0, dtype=torch.float32)
                    # for param in model.parameters():
                    #     l1_regularization += torch.norm(param, 1)
                    #     l2_regularization += torch.norm(param, 2)
                    loss =  cross_loss #+ 0.0005*l1_regularization + 0.0001*l2_regularization
                    loss.backward()
                    optimizer.step()
                    train_loss += loss
                accuracy = calculate_accuracy(model, valid_dl, False, False)
                val_loss = get_valid_loss(model, cross_loss, valid_dl)
                if val_loss < val_loss_threshold and train_loss < loss_threshold:                    
                    print("valid loss: {:.3f}, accuracy: {:.3f}".format(val_loss, accuracy))
                    flag = 0
                    break                
                
        accuracy = calculate_accuracy(model, train_dl, False, False)
        info = "train loss: {:.3f}, accuracy: {:.3f}"
        print(info.format(train_loss, accuracy))
        # predict test dataset    
        calculate_accuracy(model, test_dl, True, True)
    
    print("Predict: ", test_predicts)
    print("Label: ", test_labels)
    print("Probability: ", test_proba)
    fpr, tpr, threshold = roc_curve(test_labels, test_proba)
    auc_value = auc(fpr, tpr)
    acc = 0
    for i in range(len(test_labels)):
        if test_predicts[i] == test_labels[i]:
            acc = acc + 1
    print("AUC: {:.3f}, Accuracy: {:.3f}".format(auc_value, acc / len(test_labels)))

    

if __name__ == "__main__":
    for i in range(2): # Multiple verifications
        test_labels.clear()
        test_predicts.clear()
        test_proba.clear()
        print("the {}th: ".format(i+1))
        run()
    pass