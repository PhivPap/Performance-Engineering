import random


# CFG
universe_file = "in/in12500.tsv"
bodies = 12500
universe = ((-5e16, 5e16), (-5e16, 5e16))   # m (x, y)
mass_range = (1e10, 1e40)                   # kg
vel_range = (-1e6, 1e6)                     # m / s

# LOGIC
t = '\t'
def new_body(num):
    return {
        "id"    : str(num),
        "mass"  : str(random.uniform(mass_range[0], mass_range[1])),
        "x"     : str(random.uniform(universe[0][0], universe[0][1])),
        "y"     : str(random.uniform(universe[1][0], universe[1][1])),
        "vel_x" : str(random.uniform(vel_range[0], vel_range[1])),
        "vel_y" : str(random.uniform(vel_range[0], vel_range[1]))
    }

f = open(universe_file, "w")
f.write("id" + t + "mass" + t + "x" + t + "y" + t + "vel_x" + t + "vel_y" + "\n")

for i in range(bodies):
    b = new_body(i)
    f.write(b["id"] + t + b["mass"] + t + b["x"] + t + b["y"] + t + b["vel_x"] + t + b["vel_y"] + "\n")