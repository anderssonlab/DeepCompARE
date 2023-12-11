import torch
import torch.nn as nn

#----------------------------
# Multitasking losses
#----------------------------

def classification_loss_MT(val_pred,truth):
    """
    Method:
        Rows with only absent label (-1): ignore
        Rows with only 1 valid label: output twice the BCE loss
        Rows with 2 valid label: output the max loss of the two added by mean loss of the two
    
    Args:
        val_pred: values predicted
        truth: true binary label
        
    Output: 
        final_loss: final classification loss (scalar) 
        pred_all_correct: an indicator vector showing whether all valid labels in each row is prediced correct
    """
    # remove rows in truth with all -1
    mask=torch.logical_not(torch.all((truth==-1),dim=1))
    val_pred=val_pred[mask,:]
    truth=truth[mask,:]
    # for rows with at least one valid label
    # ignore the -1s.
    val_pred[truth==-1]=-1e9
    truth[truth==-1]=0

    loss_fun=nn.BCEWithLogitsLoss(reduction='none')
    unreduced_loss=loss_fun(val_pred,truth)
    loss_max,_=unreduced_loss.max(axis=1)
    pred_all_correct= (loss_max<0.693).int()
    final_loss=loss_max.mean()+unreduced_loss.mean()
    return final_loss,pred_all_correct



def regression_loss_MT(pred,truth,weight):
    """
    Args:
        weight: weight applied to each of the 8 tracks
    Output:
        weighted loss of each track (a vector of length 8).
    """
    
    loss_fun=nn.MSELoss(reduction="none")
    return (loss_fun(pred,truth)).mean(axis=0)*weight





#----------------------------
# Single task losses
#----------------------------

def classification_loss_ST(pred,truth):
    loss_fun=torch.nn.BCEWithLogitsLoss()
    return loss_fun(pred,truth)

def regression_loss_ST(pred,truth):
    loss_fun=torch.nn.MSELoss()
    return loss_fun(pred,truth)
