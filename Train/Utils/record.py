import torch
import pandas as pd


class Record:
    """
    A class to record regression loss (of each track) and classification accuracy (scalar).
    For regression models, self.accs=0
    For classification models, self.losses=[0]*8 or 0.
    """
   
    def __init__(self,device,n_epochs,num_tracks):
        self.losses=torch.empty((n_epochs,num_tracks),device=device) 
        self.accs=torch.empty(n_epochs,device=device)
        self.epoch=0
        self.running_data_size=0
        self.running_loss=torch.zeros(num_tracks,device=device)
        self.running_acc=0
        self.device=device
        self.num_tracks=num_tracks

    def clean(self):
        self.running_loss=torch.zeros(self.num_tracks,device=self.device)
        self.running_acc=0
        self.epoch+=1
        self.running_data_size=0


    def batch_update(self,batch_size,loss_reg=False,pred=False): 
        self.running_data_size+=batch_size
        if not isinstance(loss_reg, bool):
            self.running_loss+=loss_reg*batch_size
        if not isinstance(pred, bool):
            self.running_acc+=((pred==1).sum())

    def calculate(self):
        self.losses[self.epoch]=self.running_loss/self.running_data_size
        self.accs[self.epoch]=self.running_acc/self.running_data_size


    def return_res(self):
        self.losses=self.losses.cpu().numpy()[0:self.epoch]
        self.accs=self.accs.cpu().numpy()[0:self.epoch] 
        return self
