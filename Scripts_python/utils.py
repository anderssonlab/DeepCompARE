import numpy as np
import torch
import kipoiseq

def encode(sequence):
    return kipoiseq.transforms.functional.one_hot_dna(sequence)


def shift_seq(seq,direction,dist):
    """ 
    add "N"*dist to object seq from directin 
    """
    if direction=="right": # add 0 from left,subset from left
        return ("N"*dist+seq)[0:600]
    elif direction=="left": # add 0 from right, subset from right
        return (seq+"N"*dist)[dist:dist+600]
    else:
        raise ValueError("parameter direction not recognized!")



def find_available_gpu():
    if torch.cuda.is_available():
        # Number of GPUs available
        num_gpus = torch.cuda.device_count()

        for i in range(num_gpus):
            # Set the current GPU
            torch.cuda.set_device(i)

            # Get the name of the current GPU
            gpu_name = torch.cuda.get_device_name(i)

            # Get the current GPU's memory usage
            total_memory = torch.cuda.get_device_properties(i).total_memory
            allocated_memory = torch.cuda.memory_allocated(i)
            free_memory = total_memory - allocated_memory
            
            if free_memory>(1024**3)*20:
                return str(i)
    else:
        raise Exception("CUDA is not available.")