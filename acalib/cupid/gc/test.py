import pycupid
import numpy as np
d = np.array([1,2,3], dtype=np.double)
v = np.array([0,0,0], dtype=np.double)
cfg = {"answer": 42, "pi": np.pi, "question": "the answer to life the universe and everything"}
rms = 30
velax = 0
pycupid.gaussclumps(d,v,cfg,rms,velax)
