import torch
import numpy as np  

def regularize(angles):
    angles = (angles + torch.pi)% (2*torch.pi) - torch.pi
    return angles

def geodesic_so2(x0, x1):
    # paths = torch.stack([x1 - x0, x1 - x0 + 2*torch.pi, x1 - x0 - 2*torch.pi], dim=-1)
    # distances = paths.abs()
    # geodesics_index = torch.min(distances, dim=-1, keepdim=True).indices
    # geodesics = paths.gather(-1, geodesics_index).squeeze(-1)
    geodesics = (x1 - x0 + torch.pi) % (2*torch.pi) - torch.pi
    x_terminal = x0 + geodesics
    return geodesics, x_terminal