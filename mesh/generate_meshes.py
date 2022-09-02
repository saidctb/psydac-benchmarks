import os
cmd = "psydac-mesh --analytical quarter_annulus -n {nc} {nc} {nc} -d {d} {d} {d} -o mesh/quarter_annulus_{nc}_{d}.h5"
ncells = [64,96,128,160,192,256]
degree = [2,3,4,5]

for nc in ncells:
    for d in degree:
        os.system(cmd.format(nc=nc,d=d))

cmd = "psydac-mesh --analytical identity -n {nc} {nc} {nc} -d {d} {d} {d} -o mesh/identity_{nc}_{d}.h5"
ncells = [64,96,128,160,192,256]
degree = [2,3,4,5]

for nc in ncells:
    for d in degree:
        os.system(cmd.format(nc=nc,d=d))

