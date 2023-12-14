"""
Define models to benchmark
"""

import torch
from torch import nn
import torch.nn.functional as F
import numpy as np


def calculate_out_length(orig_length,ks,ds):
    l=orig_length
    for i in range(len(ks)):
        l=l-ds[i]*(ks[i]-1)
    return l

class Conv5D(nn.Module):
    def __init__(self,
                 ks=[15,9,9,9,9],
                 cs=[64,64,64,64,64],
                 ds=[1,2,4,8,16]
                 ):
        super().__init__()
        l=calculate_out_length(600,ks,ds)
        self.conv1=nn.Conv1d(in_channels=4, out_channels=cs[0], kernel_size=ks[0],dilation=ds[0])
        self.conv2=nn.Conv1d(in_channels=cs[0], out_channels=cs[1], kernel_size=ks[1],dilation=ds[1])
        self.conv3=nn.Conv1d(in_channels=cs[1], out_channels=cs[2], kernel_size=ks[2],dilation=ds[2])
        self.conv4=nn.Conv1d(in_channels=cs[2], out_channels=cs[3], kernel_size=ks[3],dilation=ds[3])
        self.conv5=nn.Conv1d(in_channels=cs[3], out_channels=cs[4], kernel_size=ks[4],dilation=ds[4])
        self.fc1=nn.Linear(l*cs[-1],1)
        self.d1=nn.Dropout(p=0.1)
        self.d2=nn.Dropout(p=0.1)
        self.d3=nn.Dropout(p=0.1)
        self.d4=nn.Dropout(p=0.1)
        self.d5=nn.Dropout(p=0.1)
    def forward(self,x):
        batch_size=x.size()[0]
        x=F.relu(self.conv1(x))
        x=self.d1(x)
        x=F.relu(self.conv2(x))
        x=self.d2(x)
        x=F.relu(self.conv3(x))
        x=self.d3(x)
        x=F.relu(self.conv4(x))
        x=self.d4(x)
        x=F.relu(self.conv5(x))
        x=self.d5(x)
        x=x.view(batch_size,-1)
        x=self.fc1(x)
        return x




class CRConv5D(nn.Module):
    def __init__(self,
                 ks=[15,9,9,9,9],
                 cs=[64,64,64,64,64],
                 ds=[1,2,4,8,16]
                 ):
        super().__init__()
        l=calculate_out_length(600,ks,ds)
        self.conv1=nn.Conv1d(in_channels=4, out_channels=cs[0], kernel_size=ks[0],dilation=ds[0])
        self.conv2=nn.Conv1d(in_channels=cs[0], out_channels=cs[1], kernel_size=ks[1],dilation=ds[1])
        self.conv3=nn.Conv1d(in_channels=cs[1], out_channels=cs[2], kernel_size=ks[2],dilation=ds[2])
        self.conv4=nn.Conv1d(in_channels=cs[2], out_channels=cs[3], kernel_size=ks[3],dilation=ds[3])
        self.conv5=nn.Conv1d(in_channels=cs[3], out_channels=cs[4], kernel_size=ks[4],dilation=ds[4])
        self.fc_reg=nn.Linear(l*cs[-1],1)
        self.fc_class=nn.Linear(l*cs[-1],1)
        self.d1=nn.Dropout(p=0.1)
        self.d2=nn.Dropout(p=0.1)
        self.d3=nn.Dropout(p=0.1)
        self.d4=nn.Dropout(p=0.1)
        self.d5=nn.Dropout(p=0.1)

    def forward(self,x):
        batch_size=x.size()[0]
        x=F.relu(self.conv1(x))
        x=self.d1(x)
        x=F.relu(self.conv2(x))
        x=self.d2(x)
        x=F.relu(self.conv3(x))
        x=self.d3(x)
        x=F.relu(self.conv4(x))
        x=self.d4(x)
        x=F.relu(self.conv5(x))
        rep=x.view(batch_size,-1)
        regression=self.fc_reg(rep)
        classification=self.fc_class(rep)
        return torch.cat((regression,classification),dim=1)
    
    
    


class grouped_fc_layer(nn.Module):
    def __init__(self,input_length,num_groups):
        super().__init__()
        self.weight=nn.Parameter(torch.empty(1,input_length,num_groups))
        self.bias=nn.Parameter(torch.rand(8)-0.5)
        nn.init.kaiming_normal_(self.weight)

    def forward(self,x):
        x=torch.matmul(x,self.weight)+self.bias
        x=torch.diagonal(x,dim1=1, dim2=2)
        return x










# for modeling only regression or only classification
class AstigConv5D(nn.Module):
    def __init__(self,
                 ks=[15,9,9,9,9],
                 cs=[128,128,128,128,128],
                 ds=[1,2,4,8,16],
                 gs=[1,1,1,4,8]):    
        super().__init__()
        self.l=calculate_out_length(600,ks,ds)
        self.c_group=cs[-1]//8
        self.conv1=nn.Conv1d(in_channels=4, out_channels=cs[0], kernel_size=ks[0],dilation=ds[0],groups=gs[0])
        self.conv2=nn.Conv1d(in_channels=cs[0], out_channels=cs[1], kernel_size=ks[1],dilation=ds[1],groups=gs[1])
        self.conv3=nn.Conv1d(in_channels=cs[1], out_channels=cs[2], kernel_size=ks[2],dilation=ds[2],groups=gs[2])
        self.conv4=nn.Conv1d(in_channels=cs[2], out_channels=cs[3], kernel_size=ks[3],dilation=ds[3],groups=gs[3])
        self.conv5=nn.Conv1d(in_channels=cs[3], out_channels=cs[4], kernel_size=ks[4],dilation=ds[4],groups=gs[4])
        self.d1=nn.Dropout(p=0.1)
        self.d2=nn.Dropout(p=0.1)
        self.d3=nn.Dropout(p=0.1)
        self.d4=nn.Dropout(p=0.1)
        self.d5=nn.Dropout(p=0.1)
        self.fc1=grouped_fc_layer(self.l*cs[-1]//8,8)


    def forward(self,x):
        batch_size=x.size()[0]
        x=F.relu(self.conv1(x))
        x=self.d1(x)
        x=F.relu(self.conv2(x))
        x=self.d2(x)
        x=F.relu(self.conv3(x))
        x=self.d3(x)
        x=F.relu(self.conv4(x))
        x=self.d4(x)
        x=F.relu(self.conv5(x))
        x=self.d5(x)
        x=x.reshape(batch_size,8,self.l*self.c_group)
        x=self.fc1(x)
        return x


# for modeling CR
class AstigCRConv5D(nn.Module):
    def __init__(self,
                 ks=[15,9,9,9,9],
                 cs=[128,128,128,128,128],
                 ds=[1,2,4,8,16],
                 gs=[1,1,1,4,8]):
        super().__init__()
        self.l=calculate_out_length(600,ks,ds)
        self.c_group=cs[-1]//8
        self.conv1=nn.Conv1d(in_channels=4, out_channels=cs[0], kernel_size=ks[0],dilation=ds[0],groups=gs[0])
        self.conv2=nn.Conv1d(in_channels=cs[0], out_channels=cs[1], kernel_size=ks[1],dilation=ds[1],groups=gs[1])
        self.conv3=nn.Conv1d(in_channels=cs[1], out_channels=cs[2], kernel_size=ks[2],dilation=ds[2],groups=gs[2])
        self.conv4=nn.Conv1d(in_channels=cs[2], out_channels=cs[3], kernel_size=ks[3],dilation=ds[3],groups=gs[3])
        self.conv5=nn.Conv1d(in_channels=cs[3], out_channels=cs[4], kernel_size=ks[4],dilation=ds[4],groups=gs[4])
        self.d1=nn.Dropout(p=0.1)
        self.d2=nn.Dropout(p=0.1)
        self.d3=nn.Dropout(p=0.1)
        self.d4=nn.Dropout(p=0.1)
        self.d5=nn.Dropout(p=0.1)
        self.fc_regression=grouped_fc_layer(self.l*cs[-1]//8,8)
        self.fc_classification=grouped_fc_layer(self.l*cs[-1]//8,8)

    def forward(self,x):
        batch_size=x.size()[0]
        x=F.relu(self.conv1(x))
        x=self.d1(x)
        x=F.relu(self.conv2(x))
        x=self.d2(x)
        x=F.relu(self.conv3(x))
        x=self.d3(x)
        x=F.relu(self.conv4(x))
        x=self.d4(x)
        x=F.relu(self.conv5(x))
        x=self.d5(x)
        rep=x.reshape(batch_size,8,self.l*self.c_group)
        regression=self.fc_regression(rep)
        classification=self.fc_classification(rep)
        return torch.cat((regression,classification),dim=1)






