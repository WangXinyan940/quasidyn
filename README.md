# quasidyn

This program is for quasi-classical direct dynamics. The conformation is chosen as transition state. Each positive normal mode of the TS is given its zero point energy plus a random energy from Boltzmann distribution. The sign of velocity is given randomly. For negative vibration mode, a random energy from Blotzmann sampling is given. The sign is setted in configure file and can be "forward", "reverse" and "both". The direction of virtual vibration mode can be seen by the tool "showdirection.py". The trajectory is integrated by Velocity Verlet integration. 

## Required

- QM engine
  - Gaussian
  - ORCA (comming soon)
  
## Usage

    python quasidyn.py [configure file]

If not setting the name of configure file, "config.json" will be used by default.

## Config file example:

    {
            "QM_engine":"GAUSSIAN", #Only Gaussian is supported at present.
            "nproc":"10", #Number of processors.
            "mem":"2000mw", #Memory used.
            "method":"B3LYP em=gd3", #Computational method.
            "basis":"6-31G* int=fine", #Basis set.
            "extra":"", #Extra setting needed, such as "scf" setting.
            "charge":-1, #System charge.
            "multi":1, #Spin multiplicity.
            "call_gaussian":"g16", #The command needed to call Gaussian.

            "initfile":"ts.fchk", #The name of init file. A fchk file from frequency analysis is needed. 
                                  #Fchk file can be transferred from chk file by formchk tool from Gaussian.
            "temperature":300, #Temperature (K).
            "direction":"forward", #Direction of virtual vibration mode. [forward/reverse/both]
            "traj_filename":"traj.xyz", #Name of trajectory file. The energy of each step is also written in trajectory file.
            "dt":0.5, #Timestep (fs). The timestep is always 0.5 fs. For big organic molecule, 1.0 fs is sometimes used.
            "tmax":200, #Max length of trajectory (fs). 
            "maxe":5, #The maximum energy change between two steps (kcal/mol). If energy change is larger than this value,
                      #the step need to be seen as unphysical and need to stop. 
            "states":{
            #Describe states on potential energy surface. If the system goes to one of these states, 
            #it means the trajectory touchs the end and need to be stopped.
                "state-A":[
            {
                "type":"A", #Can be "A" (angle) and "B" (bond). "A" is the angle between 
                            #vector ab and vector ac. "B" is the length of vector ab.
                "ia":1, #The index of atom a. Start from 1.
                "ib":4,
                "ic":5,
                "min":0.0, #The minimum value. The unit is radian for "A" and angstrom for "B".
                "max":2.2
            }
        ],
        "state-B":[
            {
                "type":"B",
                "ia":1,
                "ib":6,
                "min":0.0,
                "max":2.2
            }
        ]
            },
            "writevel":false #If velocity is written. The unit in velocity file is angstrom/fs.
    }
