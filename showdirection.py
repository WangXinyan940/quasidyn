import quasidyn
import json
import numpy as np

with open("config.json", "r") as f:
    config = json.load(f)

e, V = quasidyn.readGaussVib(config["initfile"])
atoms, xyz, mass = quasidyn.readGaussBasic(config["initfile"])
for n,ei in enumerate(np.argsort(np.abs(e))):
    if ei < 6:
        continue
    if e[n] < 0:
        break
print("Virtual Freq: {} cm^-1".format(np.sqrt(np.abs(e[n])) / 2. / np.pi / 2.9979e10))
vxyz = np.matrix(xyz.reshape((xyz.shape[0] * xyz.shape[1], 1)))
for ni in [0.0,0.1,0.2,0.3,0.4,0.5]:
    ef = vxyz + V[:,n] * ni * quasidyn.ANGSTROM
    er = vxyz - V[:,n] * ni * quasidyn.ANGSTROM
    xyz_f = ef.reshape(xyz.shape)
    xyz_r = er.reshape(xyz.shape)

    quasidyn.writeTraj(0, atoms, xyz_f, 0, None, "c1.xyz", writevel=False)
    quasidyn.writeTraj(0, atoms, xyz_r, 0, None, "c-1.xyz", writevel=False)