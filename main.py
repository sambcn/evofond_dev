import fractions
import numpy as np
import matplotlib.pyplot as plt
from time import time
from src.utils import hydrogrammeLavabre, read_hecras_data, reverse_data, add_injections
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
from src.sedimentTransport.meunier1989 import Meunier1989
from src.sedimentTransport.meyerpeter1948 import MeyerPeter1948

def main():

    ############################################
    # granulometry = Granulometry(0.07, 0.02, 0.05, 0.2, 0.08, 0.15, 2.6)
    # profile = build_profile([0, 200, 400, 600, 800], [195, 395, 595, 795, 1000], 5, [6, 8, 6, 8, 7], [0.1, 0.05, 0.001, 0.07, 0.001], granulometry, z0=400, height_drops=[3, 0, -4, 0, 2], z_min_diff=2, manning=None, section="Trapezoidal")
    ############################################
    granulometry = Granulometry(0.07, 0.02, 0.05, 0.2, 0.08, 0.15, 2.6)
    profile = build_profile([0, 200], [190, 400], 5, [8, 4], [0.005, 0.05], granulometry, z0=400, height_drops=[0, 0], z_min_diff=3, manning=None, put_exit_loss=True, section="Trapezoidal")
    ############################################ 
    # granulometry = Granulometry(0.07, 0.02, 0.05, 0.2, 0.08, 0.15, 2.6)
    # profile = build_profile([0, 200], [190, 390], 10, [8, 4], [0.05, 0.05], granulometry, z0=400, height_drops=[0, 0], z_min_diff=1, manning=None, put_exit_loss=True, section="Trapezoidal", b0_min=[8, 4])
    # profile2 = build_profile([0, 200], [190, 390], 10, [8, 4], [0.05, 0.05], granulometry, z0=400, height_drops=[0, 0], z_min_diff=1, manning=None, put_exit_loss=True, section="Rectangular", b0_min=[8, 4])
    ############################################
    # granulometry = Granulometry(0.07, 0.02, 0.05, 0.2, 0.08, 0.15, 2.6)
    # profile = build_profile([0, 200], [190, 390], 10, [8, 8], [0.05, 0.001], granulometry, z0=400, height_drops=[0, 0], z_min_diff=1, manning=None, put_exit_loss=True, section="Rectangular")
    ############################################
    
    profile.complete(10)

    Q_list = [15 for i in range(profile.get_nb_section())]
    injections = [0 + (10 if i==20 else 0) + (-5 if i == 50 else 0) for i in range(profile.get_nb_section())]
    Q_list = add_injections(Q_list, injections)
    y_list = profile.compute_depth(Q_list=Q_list, plot=True, friction_law="Manning-Strickler", upstream_condition="critical_depth", downstream_condition="critical_depth", method="ImprovedEuler")

    # t = np.arange(0, 3600, 10)
    # hydrogram = hydrogrammeLavabre(30, 1200, 3, 10, t)
    # sedimentogram = hydrogram*0
    # law = Rickenmann1990()
    # friction_law = "Ferguson"
    # frac_of_simu = 1
    # cfl=0.99
    # critical=True
    # upstream_condition="critical_depth"
    # downstream_condition="critical_depth"
    # plot=True
    # animate=False
    # whatever = profile.compute_event(hydrogram[:int(len(hydrogram)*frac_of_simu)], t[:int(len(hydrogram)*frac_of_simu)], law, sedimentogram=sedimentogram, friction_law=friction_law, cfl=cfl, critical=critical, upstream_condition=upstream_condition, downstream_condition=downstream_condition, plot=plot, animate=animate)

    return



def build_profile(x_begin_list, x_end_list, dx, widths, slopes, granulometries, z0=200, height_drops=None, z_min_diff=0, manning=None, tauc_over_rho=1, K_over_tauc=0.3, section="Rectangular", put_exit_loss=True, b0_min=None):
    """
    Made by and for devs.
    """
    if type(granulometries) == Granulometry:
        granulometries = [granulometries for _ in range(len(widths))]
    x = []
    limit_index = []
    for i in range(len(x_begin_list)):
        x += list(np.arange(x_begin_list[i], x_end_list[i]-0.5*dx, dx)) + [x_end_list[i]]
        limit_index.append(len(x))
    n = len(x)
    height_drops = [0 for _ in range(len(x_begin_list)-1)] if height_drops==None else height_drops
    z = [z0]
    b = []
    b_min = []
    granulometry = []
    exit_loss_coef_list = []
    old_part_index = 0
    for i in range(n):
        part_index = old_part_index if i < limit_index[old_part_index] else old_part_index+1
        b.append(widths[part_index])
        if section=="Trapezoidal":
            b_min.append(b0_min[part_index] if b0_min != None else b[i])
        if i > 0 and b[-1] != b[-2] and put_exit_loss:
            exit_loss_coef_list.append(1)
        else:
            exit_loss_coef_list.append(0)
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
            s = RectangularSection(x[i], z[i], b[i], y_max=10, z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        elif section == "Trapezoidal":
            s = TrapezoidalSection(x[i], z[i], b[i], b_min[i], f_left=0.3, f_right=0.5, y_max=4, z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        elif section == "Irregular":
            points = [(0, 50), (0, 0), (b[i], 0), (b[i], 50)] # rectangular, but it will be treated as an irregular section
            s = IrregularSection(points, x[i], z[i], z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        section_list.append(s)
    profile = Profile(section_list, exit_loss_coef_list=exit_loss_coef_list)
    return profile























if __name__=="__main__":
    Performance.start()
    ani = main()
    Performance.stop()
    Performance.print_perf()
    # Performance.save_perf("./results/perf/perf.txt")
    plt.show()