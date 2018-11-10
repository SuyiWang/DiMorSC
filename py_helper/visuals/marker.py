from vispy import visuals, scene, color
import numpy as np

def new(tc, pos):
    mkr = scene.visuals.Markers(parent=tc)
    mkr.set_data(pos)
    return mkr

def update(mkr, pos):
    mkr.set_data(pos = p)
