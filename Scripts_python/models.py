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









#---------------------------------
# BPNet
#----------------------------------


class ScalarHead21(nn.Module):
    def __init__(self):        
        super().__init__()
        self.global_avg_pool = nn.AdaptiveAvgPool1d(1)
        # Define hidden layers
        self.hidden_layers = nn.Sequential(
            nn.Linear(21, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU(),
            nn.Dropout(0.1)
        )
        
        self.final_layer = nn.Linear(32, 1)

    def forward(self, x):
        x = self.global_avg_pool(x)
        x = x.squeeze() # [B, C]
        x = self.hidden_layers(x)
        x = self.final_layer(x)
        return x


class DilatedConv1D21(nn.Module):
    def __init__(
            self,
            conv1_kernel_size=25,
            n_dil_layers=6,
    ):
        super().__init__()
        self.conv1_kernel_size = conv1_kernel_size
        self.n_dil_layers = n_dil_layers
        self.first_conv = nn.Conv1d(4, 21, 25, padding='same') 
        self.dil_conv_layers = nn.ModuleList()
        for i in range(1, self.n_dil_layers + 1):
            dil_conv = nn.Conv1d(21,21,3, padding='same', dilation=2**i) # no compression along the layers
            self.dil_conv_layers.append(dil_conv)

    def forward(self, x):
        x = F.relu(self.first_conv(x))
        prev_layer = x
        for dil_conv in self.dil_conv_layers:
            x = prev_layer
            conv_output = F.relu(dil_conv(x))
            prev_layer = prev_layer + conv_output # residual connect
        return prev_layer
    

class BPNet21(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = DilatedConv1D21()
        self.head = ScalarHead21()

    def forward(self, x):
        x = self.net(x)
        x = self.head(x)
        return x
    
    
    
    
    
    
    
    
    


    
    

class ScalarHead64(nn.Module):
    def __init__(self):        
        super().__init__()
        self.global_avg_pool = nn.AdaptiveAvgPool1d(1)
        # Define hidden layers
        self.hidden_layers = nn.Sequential(
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU(),
            nn.Dropout(0.1)
        )
        
        self.final_layer = nn.Linear(32, 1)

    def forward(self, x):
        x = self.global_avg_pool(x)
        x = x.squeeze() # [B, C]
        x = self.hidden_layers(x)
        x = self.final_layer(x)
        return x


class DilatedConv1D64(nn.Module):
    def __init__(
            self,
            conv1_kernel_size=25,
            n_dil_layers=6,
    ):
        super().__init__()
        self.conv1_kernel_size = conv1_kernel_size
        self.n_dil_layers = n_dil_layers
        self.first_conv = nn.Conv1d(4, 64, 25, padding='same') 
        self.dil_conv_layers = nn.ModuleList()
        for i in range(1, self.n_dil_layers + 1):
            dil_conv = nn.Conv1d(64,64,3, padding='same', dilation=2**i) # no compression along the layers
            self.dil_conv_layers.append(dil_conv)

    def forward(self, x):
        x = F.relu(self.first_conv(x))
        prev_layer = x
        for dil_conv in self.dil_conv_layers:
            x = prev_layer
            conv_output = F.relu(dil_conv(x))
            prev_layer = prev_layer + conv_output # residual connect
        return prev_layer
    

class BPNet64(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = DilatedConv1D64()
        self.head = ScalarHead64()

    def forward(self, x):
        x = self.net(x)
        x = self.head(x)
        return x




