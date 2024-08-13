import torch
from loss import *
from record import Record
import math
from abc import ABC, abstractmethod

def load_data_and_predict(model, dat, idx, device):
    """ 
    Loads npyDataset and performs model prediction.
    Only used for training and validation
    """
    x = torch.tensor(dat[idx][0], device=device)
    return model(x)


# Run model and calculate loss for different tasks
# Output order:
# loss_backprop, loss_reg_to_record (a vector of length 8), pred_classification_to_record
def run_CR_MT(device, model, dat, idx, weight):
    out = load_data_and_predict(model, dat, idx, device)
    y = torch.tensor(dat.Y[idx], device=device)
    z = torch.tensor(dat.Z[idx], device=device)
    loss_reg = regression_loss_MT(out[:, 0:8], z, weight)
    loss_class, pred_all_correct = classification_loss_MT(out[:, 8:16], y)
    loss_backprop=loss_reg.mean()+loss_class
    return loss_backprop, loss_reg.detach(), pred_all_correct

def run_regression_MT(device, model, dat, idx, weight):
    out = load_data_and_predict(model, dat, idx, device)
    z = torch.tensor(dat.Z[idx], device=device)
    loss_reg = regression_loss_MT(out, z, weight)
    return loss_reg.mean(),loss_reg.detach(),False

def run_classification_MT(device, model, dat, idx, _):
    """
    Args:
        _ is a place holder to make sure all run_XXX models take in 5 params
    """
    out = load_data_and_predict(model, dat, idx, device)
    y = torch.tensor(dat.Y[idx], device=device)
    loss_class, pred_all_correct=classification_loss_MT(out, y)
    return loss_class, False, pred_all_correct





# Output order:
# loss_backprop, loss_reg_to_record (a vector of length 8), pred_classification_to_record

def run_CR_ST(device, model, dat, idx, _):
    out = load_data_and_predict(model, dat, idx, device)
    y = torch.tensor(dat[idx][1], device=device)
    z = torch.tensor(dat[idx][2], device=device)
    loss_reg=regression_loss_ST(out[:,0:1], z)
    loss_class=classification_loss_ST(out[:, 1:2], y)
    loss_backprop=loss_reg+loss_class
    pred=(out[:, 1:2]>0)
    return loss_backprop, loss_reg.detach(), (pred==y).long()

def run_regression_ST(device, model, dat, idx, _):
    out = load_data_and_predict(model, dat, idx, device)
    z = torch.tensor(dat[idx][2], device=device)
    loss_reg = regression_loss_ST(out, z)
    return loss_reg,loss_reg.detach(),False

def run_classification_ST(device, model, dat, idx, _):
    """
    Args:
        _ is a place holder to make sure all run_XXX models take in 5 params
    """
    out = load_data_and_predict(model, dat, idx, device)
    y = torch.tensor(dat[idx][1], device=device)
    loss_class=classification_loss_ST(out, y)
    pred=(out>0).long()
    return loss_class, False, (pred==y).long()





#-------------------------
# Trainer framework
#-------------------------

class Trainer(ABC):
    """ 
    A abstract class for Trainer.
    """
    def __init__(self, device, model, train_dat, val_dat, param_dict, optimizer):
        self.device = device
        self.model = model
        self.train_dat = train_dat
        self.val_dat = val_dat
        self.optimizer = optimizer
        self.n_epochs = param_dict["n_epochs"]
        self.batch_size = param_dict["batch_size"]
        self.batch_size_val = self.batch_size * 2


    def _train_epoch(self, train_function):
        self.model.train()
        idx_array = torch.randperm(self.train_dat.len, device=self.device)
        for i in range(math.ceil(self.train_dat.len / self.batch_size)):
            start_idx, end_idx = self.batch_size*i, min(self.batch_size*(i+1), self.train_dat.len)
            idx = idx_array[start_idx:end_idx].cpu().numpy()  
            self.optimizer.zero_grad()
            loss_backprop, loss_reg, pred_all_correct = train_function(self.device, self.model, self.train_dat, idx, self.weight)
            loss_backprop.backward()
            self.optimizer.step()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=5, norm_type=2.0)
            self.record_train.batch_update(len(idx),loss_reg,pred_all_correct)

        self.record_train.calculate()
        self.record_train.clean()

    def _val(self, val_function):
        self.model.eval()
        with torch.no_grad():
            for i in range(math.ceil(self.val_dat.len / self.batch_size_val)):
                start_idx, end_idx = self.batch_size_val*i, min(self.batch_size_val*(i+1), self.val_dat.len)
                idx = list(range(start_idx, end_idx))
                _, loss_reg, pred_all_correct = val_function(self.device, self.model, self.val_dat, idx, self.weight)
                self.record_val.batch_update(len(idx),loss_reg,pred_all_correct)

            self.record_val.calculate()
            self.record_val.clean()
    
    def train(self):
        for epoch in range(self.n_epochs):
            self._train_epoch()
            if epoch % 2 == 0:
                self.val()

        return self.record_train.return_res(), self.record_val.return_res()
    




#-------------------------
# Multitasking trainers
#-------------------------

class CR_MT(Trainer):
    """
    Trainer for classification+regression tasks, 8 tracks
    """
    def __init__(self, device, model, train_dat, val_dat, param_dict, optimizer):
        super().__init__(device, model, train_dat, val_dat, param_dict, optimizer)
        self.weight = torch.tensor(param_dict["weight"], device=device)
        self.record_train = Record(device, self.n_epochs,8)
        self.record_val = Record(device, self.n_epochs,8)

    def train_epoch(self):
        self._train_epoch(run_CR_MT)

    def val(self):
        self._val(run_CR_MT)

class Reg_MT(Trainer):
    """
    Trainer for regression tasks, 8 tracks
    """
    def __init__(self, device, model, train_dat, val_dat, param_dict, optimizer):
        super().__init__(device, model, train_dat, val_dat, param_dict, optimizer)
        self.weight = torch.tensor(param_dict["weight"], device=device)
        self.record_train = Record(device, self.n_epochs,8)
        self.record_val = Record(device, self.n_epochs,8)

    def train_epoch(self):
        self._train_epoch(run_regression_MT)

    def val(self):
        self._val(run_regression_MT)


class Class_MT(Trainer):
    """
    Trainer for classification tasks, 8 tracks
    """
    def __init__(self, device, model, train_dat, val_dat, param_dict, optimizer):
        super().__init__(device, model, train_dat, val_dat, param_dict, optimizer)
        self.weight = False
        self.record_train = Record(device, self.n_epochs,8)
        self.record_val = Record(device, self.n_epochs,8)

    def train_epoch(self):
        self._train_epoch(run_classification_MT)

    def val(self):
        self._val(run_classification_MT)






#-------------------------
# Singletasking trainers
#-------------------------
class CR_ST(Trainer):
    """
    Trainer for classification+regression tasks, 8 tracks
    """
    def __init__(self, device, model, train_dat, val_dat, param_dict, optimizer):
        super().__init__(device, model, train_dat, val_dat, param_dict, optimizer)
        self.weight = False
        self.record_train = Record(device, self.n_epochs,1)
        self.record_val = Record(device, self.n_epochs,1)

    def train_epoch(self):
        self._train_epoch(run_CR_ST)

    def val(self):
        self._val(run_CR_ST)


class Reg_ST(Trainer):
    """
    Trainer for regression tasks, 8 tracks
    """
    def __init__(self, device, model, train_dat, val_dat, param_dict, optimizer):
        super().__init__(device, model, train_dat, val_dat, param_dict, optimizer)
        self.weight = False
        self.record_train = Record(device, self.n_epochs,1)
        self.record_val = Record(device, self.n_epochs,1)

    def train_epoch(self):
        self._train_epoch(run_regression_ST)

    def val(self):
        self._val(run_regression_ST)


class Class_ST(Trainer):
    """
    Trainer for classification tasks, 8 tracks
    """
    def __init__(self, device, model, train_dat, val_dat, param_dict, optimizer):
        super().__init__(device, model, train_dat, val_dat, param_dict, optimizer)
        self.weight = False
        self.record_train = Record(device, self.n_epochs,1)
        self.record_val = Record(device, self.n_epochs,1)

    def train_epoch(self):
        self._train_epoch(run_classification_ST)

    def val(self):
        self._val(run_classification_ST)

