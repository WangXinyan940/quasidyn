import numpy as np 
import sys, os
import json

BOHR = 5.291772108e-11#Bohr -> m
ANGSTROM = 1e-10#angstrom -> m
AMU = 1.660539040e-27#amu -> kg
FS = 1e-15#fs -> s
EH = 4.35974417e-18#Hatree -> J
H = 6.626069934e-34
KB = 1.38064852e-23
ATOM_TO_MASS = {1:1.008, "H":1.008,
                6:12.011, "C":12.011,
                7:14.007, "N":14.007,
                8:15.999, "O":15.999,
                17:35.453, "Cl":35.453,
                35:79.904, "Br":79.904,}
#read config
with open("config.json", "r") as f:
    config = json.load(f)

def findline(word, text):
    for n,i in enumerate(text):
        if word in i:
            return n
    return -1

def readGaussBasic(fname):
    #return atoms, xyz(A), mass(A.U.)
    with open(fname, "r") as f:
        text = f.readlines()

    p = findline("Atomic numbers", text)
    totalnum = int(text[p].strip().split()[-1])
    nlines = totalnum // 6
    if totalnum % 6 > 0:
        nlines += 1
    tmp = [[int(j) for j in i.strip().split()] for i in text[p+1:p+1+nlines]]
    atoms = []
    for i in tmp:
        for j in i:
            atoms.append(j)

    p = findline("Current cartesian coordinates", text)
    totalnum = int(text[p].strip().split()[-1])
    nlines = totalnum // 5
    if totalnum % 5 > 0:
        nlines += 1
    tmp = [[float(j) for j in i.strip().split()] for i in text[p+1:p+1+nlines]]
    xyz = []
    for i in tmp:
        for j in i:
            xyz.append(j)
    xyz = np.array(xyz).reshape((len(atoms), 3)) * BOHR
    mass = np.array([ATOM_TO_MASS[i] * AMU for i in atoms])
    return atoms, xyz, mass

def readGaussGrad(fname,natoms):
    #return energy(Hatree), gradient(Hatree/Bohr)
    with open("force.log", "r") as f:
        text = f.readlines()
        ener = [i for i in text if "SCF Done:" in i]
        if len(ener) != 0:
            ener = ener[-1]
            ener = np.float64(ener.split()[4])
        else:
           ener = np.float64([i for i in text if "Energy=" in i][-1].split()[1]) 
        for ni, li in enumerate(text):
            if "Forces (Hartrees/Bohr)" in li:
                break
        forces = text[ni+3:ni+3+natoms]
        forces = [i.strip().split()[-3:] for i in forces]
        forces = [[np.float64(i[0]), np.float64(i[1]), np.float64(i[2])] for i in forces]
    return ener * EH,-np.array(forces) * EH / BOHR

def readGaussVib(fname):
    #read Hess from fchk file
    with open(fname, "r") as f:
        text = f.readlines()

    p = findline("Atomic numbers", text)
    totalnum = int(text[p].strip().split()[-1])
    nlines = totalnum // 6
    if totalnum % 6 > 0:
        nlines += 1
    tmp = [[int(j) for j in i.strip().split()] for i in text[p+1:p+1+nlines]]
    atoms = []
    for i in tmp:
        for j in i:
            atoms.append(j)

    M = np.zeros((3 * len(atoms), 3 * len(atoms)))
    M = np.matrix(M,dtype=np.float64)

    for n,a in enumerate(atoms):
        M[3*n,3*n] = ATOM_TO_MASS[a] ** -0.5
        M[3*n + 1,3*n + 1] = ATOM_TO_MASS[a] ** -0.5
        M[3*n + 2,3*n + 2] = ATOM_TO_MASS[a] ** -0.5

    p = findline("Cartesian Force Constants", text)
    totalnum = int(text[p].strip().split()[-1])
    nlines = totalnum // 5
    if totalnum % 5 > 0:
        nlines += 1

    hess = [[np.float64(j) for j in i.strip().split()] for i in text[p+1:p+1+nlines]]
    H = np.zeros((3 * len(atoms), 3 * len(atoms)))
    H = np.matrix(H,dtype=np.float64)

    pi, pj = 0, 0
    for line in hess:
        for item in line:
            H[pi,pj] = item
            H[pj,pi] = item
            pj += 1
            if pj > pi:
                pi += 1
                pj = 0

    F = M * H * M
    F = F * EH / (BOHR ** 2 * AMU)
    e,V = np.linalg.eig(F)
    return e,V

def genGaussGrad(atoms):
    def calcGrad(xyz):
        nxyz = xyz / ANGSTROM
        with open("force.gjf", "w") as f:
            f.write("%nproc={}\n".format(config["nproc"]))
            f.write("%mem={}\n".format(config["mem"]))
            f.write("%chk=force.chk\n")
            if "force.chk" in os.listdir("."):
                f.write("#p {} {} nosymm force guess=read {}\n\n".format(config["method"], config["basis"], config["extra"]))
            else:
                f.write("#p {}/{} nosymm force {}\n\n".format(config["method"], config["basis"], config["extra"]))
            f.write("Title\n\n{} {}\n".format(config["charge"], config["multi"]))
            for ni in range(len(atoms)):
                f.write("{} {} {} {}\n".format(atoms[ni], nxyz[ni,0], nxyz[ni,1], nxyz[ni,2]))
            f.write("\n\n\n\n\n")
        os.system("{} force.gjf".format(config["call_gaussian"]))
        energy, grad = readGaussGrad("force.log",len(atoms))
        return energy, grad
    return calcGrad

def vv_step(xyz, vel, grad, calcGrad, mass, dt):
    xyz_1 = xyz + vel * (dt * FS) - grad / 2.0 / mass.reshape((mass.shape[0],1)) * (dt * FS) ** 2
    energy_1, grad_1 = calcGrad(xyz_1)
    vel_1 = vel - (grad_1 + grad) / 2.0 / mass.reshape((mass.shape[0],1)) * (dt * FS)
    return xyz_1, vel_1, energy_1, grad_1

def getBoltzmann(temperature=300): #return J
    return - 0.5 * KB * temperature * np.log(np.random.random()) 

def initByVib(xyz, mass, force_consts, vibs, temperature=300, virtsign=0):
    nxyz = xyz + (np.random.random(xyz.shape) - 0.5) * 2 * 0.02 * ANGSTROM
    freq = np.sign(force_consts) * np.sqrt(np.abs(force_consts)) / 2. / np.pi
    #return velf, velr
    mv = np.zeros((len(mass) * 3,))
    for n,i in enumerate(mass):
        mv[3*n:3*n+3] = i
    velf, velr = np.zeros(xyz.shape), np.zeros(xyz.shape)
    fc = [[n,freq[n]] for n,i in enumerate(np.argsort(np.abs(freq))) if i > 5]
    for n,f in fc:
        if f < 0:
            energy = getBoltzmann(temperature=temperature)
            alpha = np.sqrt(2 * energy / (mv * np.power(vibs[:,n], 2)).sum())
            if virtsign == 0:
                velf = velf + (alpha * vibs[:,n]).reshape(xyz.shape)
                velr = velr - (alpha * vibs[:,n]).reshape(xyz.shape)
            else:
                velf = velf + virtsign * (alpha * vibs[:,n]).reshape(xyz.shape)
                velr = velr + virtsign * (alpha * vibs[:,n]).reshape(xyz.shape)
            #v = alpha * vibs[:,n]
            #print(f/2.9979e10,(mv * np.power(v,2)).sum() / KB)
        else:
            energy = 0.5 * f * H + getBoltzmann(temperature=temperature) #ZPE + K
            #energy = getBoltzmann(temperature=temperature)
            #print(energy)
            sign = -1 if np.random.randint(0,2) == 0 else 1
            alpha = np.sqrt(2 * energy / (mv * np.power(vibs[:,n], 2)).sum())
            velf = velf + (sign * alpha * vibs[:,n]).reshape(xyz.shape)
            velr = velr + (sign * alpha * vibs[:,n]).reshape(xyz.shape)
            #v = alpha * vibs[:,n]
            #print(f/2.9979e10,(mv * np.power(v,2)).sum() / KB)
    return nxyz, velf, velr

def initCartesian(xyz, mass, temperature=300):
    vel = np.zeros(xyz.shape)
    for i in range(xyz.shape[0]):
        for j in range(xyz.shape[1]):
            energy = getBoltzmann(temperature=temperature)
            sign = -1 if np.random.randint(0,2) == 0 else 1
            vel[i,j] = sign * np.sqrt(2 * energy / mass[i])
    return xyz, vel

def writeTraj(step, atoms, xyz, energy, vel, fname, writevel=False):
    with open(fname, "a") as f:
        f.write("%i\n"%xyz.shape[0])
        f.write("STEP%i: %f\n"%(step, energy / EH))
        for l in range(len(atoms)):
            f.write("%s%12.8f%12.8f%12.8f\n"%(atoms[l], xyz[l,0] / ANGSTROM, xyz[l,1] / ANGSTROM, xyz[l,2] / ANGSTROM))
    if writevel:
        with open("vel-"+fname, "a") as f:
            f.write("%i\n"%xyz.shape[0])
            f.write("STEP%i: %f\n"%(step, energy))
            for l in range(len(atoms)):
                f.write("%s%12.8f%12.8f%12.8f\n"%(atoms[l], vel[l,0] / (ANGSTROM / FS), vel[l,1] / (ANGSTROM / FS), vel[l,2] / (ANGSTROM / FS)))

def ifStop(xyz, states):
    for state in states.values():
        if len(state) == 0:
            continue
        ifstate = []
        for cond in state:
            if cond["type"] == "B":
                ia,ib,dmin,dmax = cond["ia"], cond["ib"], cond["min"], cond["max"]
                dist = np.sqrt(np.power(xyz[ia-1,:] - xyz[ib-1,:],2).sum()) / ANGSTROM
                if dist < dmax and dist > dmin:
                    ifstate.append(True)
                else:
                    ifstate.append(False)
            if cond["type"] == "A":
                ia,ib,ic,amin,amax = cond["ia"], cond["ib"], cond["ic"], cond["min"], cond["max"]
                vi = xyz[ia-1,:] - xyz[ib-1,:]
                vj = xyz[ia-1,:] - xyz[ic-1,:]
                angle = np.argcos(np.dot(vi,vj) / np.sqrt((vi ** 2).sum()) / np.sqrt((vj ** 2).sum()))
                if angle < amax and angle > amin:
                    ifstate.append(True)
                else:
                    ifstate.append(False)
        if False not in ifstate:
            return True
    return False

if __name__ == '__main__':
    if config["QM_engine"] == "GAUSSIAN":
        readBasic = readGaussBasic
        genGrad = genGaussGrad
        readGrad = readGaussGrad
        readVib = readGaussVib

    atoms, xyz, mass = readBasic(config["initfile"])
    calcGrad = genGrad(atoms)
    force_consts, vibs = readVib(config["initfile"])
    if config["direction"] == "forward":
        xyz, vel, velr = initByVib(xyz, mass, force_consts, vibs, temperature=config["temperature"], virtsign=1)
    elif config["direction"] == "reverse":
        xyz, vel, velr = initByVib(xyz, mass, force_consts, vibs, temperature=config["temperature"], virtsign=-1)
    elif config["direction"] == "both":
        xyz, vel, velr = initByVib(xyz, mass, force_consts, vibs, temperature=config["temperature"], virtsign=0)
    energy, grad = calcGrad(xyz)
    if config["direction"] == "both":
        init_xyz = np.zeros(xyz.shape)
        init_vel = np.zeros(vel.shape)
        init_grad = np.zeros(grad.shape)
        init_xyz[:,:] = xyz[:,:]
        init_vel[:,:] = velr[:,:]
        init_grad[:,:] = grad[:,:]
        init_energy = float(energy)
    step = 0
    energy_list = []
    while True:
        step += 1
        if ifStop(xyz, config["states"]):
            print("ARRIVE AT ONE STATE. STOP.")
            break
        if step * config["dt"] > config["tmax"]:
            break
        writeTraj(step, atoms, xyz, energy, vel, config["traj_filename"], writevel=config["writevel"])
        energy_list.append(energy)
        xyz, vel, energy, grad = vv_step(xyz, vel, grad, calcGrad, mass, config["dt"])
        if np.abs(energy_list[-1] - energy) * 627.5 > config["maxe"]:
            print("ELECTRON ENERGY CHANGE TOO MUCH!!!")
            break
    if not config["direction"] == "both":
        exit()
    step = 0
    energy_list = []
    xyz = init_xyz
    vel = init_vel
    grad = init_grad
    energy = init_energy
    while True:
        step += 1
        if ifStop(xyz, config["states"]):
            print("ARRIVE AT ONE STATE. STOP.")
            break
        if step * config["dt"] > config["tmax"]:
            break
        writeTraj(step, atoms, xyz, energy, vel, "R-"+config["traj_filename"], writevel=config["writevel"])
        energy_list.append(energy)
        xyz, vel, energy, grad = vv_step(xyz, vel, grad, calcGrad, mass, config["dt"])
        if np.abs(energy_list[-1] - energy) * 627.5 > config["maxe"]:
            print("ELECTRON ENERGY CHANGE TOO MUCH!!!")
            break