# quasidyn

This program is for quasi-classical direct dynamics. The conformation is chosen as transition state. Each positive normal mode of the TS is given its zero point energy plus a random energy from Boltzmann distribution. The sign of velocity is given randomly. For negative vibration mode, a random energy from Blotzmann sampling is given. The sign is setted in configure file and can be "forward", "reverse" and "both". The direction of virtual vibration mode can be seen by the tool "showdirection.py".

Usage:
    python quasidyn.py [config name]

If not set configure file name, "config.json" will be read by default.

