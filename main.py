import numpy as np
import matplotlib.pyplot as plt
from time import time
from src.utils import hydrogrammeLavabre, read_hecras_data, reverse_data
from src.irregularSection import IrregularSection
from src.rectangularSection import RectangularSection
from src.trapezoidalSection import TrapezoidalSection
from src.profile import Profile
from src.granulometry import Granulometry
from src.perf import Performance
from src.sedimentTransport.sedimentTransportLaw import SedimentTransportLaw 
from src.sedimentTransport.rickenmann1990 import Rickenmann1990 
from src.sedimentTransport.rickenmann1991 import Rickenmann1991
from src.sedimentTransport.lefort2015 import Lefort2015
from src.sedimentTransport.lefortsogreah1991 import LefortSogreah1991


def main():
    print("Helloworld")
    return



def build_profile(x_begin_list, x_end_list, dx, widths, slopes, granulometries, z0=200, height_drops=None, z_min_diff=0, manning=0.013, tauc_over_rho=1, K_over_tauc=0.3, section="Rectangular"):
    """
    Made by and for devs.
    """
    x = []
    limit_index = []
    for i in range(len(x_begin_list)):
        x += list(np.arange(x_begin_list[i], x_end_list[i]-0.5*dx, dx)) + [x_end_list[i]]
        limit_index.append(len(x))
    n = len(x)
    height_drops = [0 for _ in range(len(x_begin_list)-1)] if height_drops==None else height_drops
    z = [z0]
    b = []
    granulometry = []
    old_part_index = 0
    for i in range(n):
        part_index = old_part_index if i < limit_index[old_part_index] else old_part_index+1
        b.append(widths[part_index])
        granulometry.append(granulometries[part_index])
        if i == 0:
            old_part_index = part_index
            continue
        else:
            z_new = z[-1]-(x[i]-x[i-1])*slopes[old_part_index]
            if part_index > old_part_index:
                z_new -= height_drops[old_part_index]
            z.append(z_new)
            old_part_index = part_index

    section_list = []
    for i in range(n):
        if section == "Rectangular":
            s = RectangularSection(x[i], z[i], b[i], y_max=50, z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        elif section == "Trapezoidal":
            s = TrapezoidalSection(x[i], z[i], b[i], 1, y_max=10, z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        elif section == "Irregular":
            points = [(0, 50), (0, 0), (b[i], 0), (b[i], 50)] # rectangular, but it will be treated as an irregular section
            s = IrregularSection(points, x[i], z[i], z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        section_list.append(s)
    profile = Profile(section_list)
    return profile























if __name__=="__main__":
    # Performance.start()
    ani = main()
    # Performance.stop()
    # Performance.print_perf()
    # Performance.save_perf("./results/perf/perf.txt")
    plt.show()