from silo_calcs import SiloCalcs
import plot_graphs

opc1 = SiloCalcs(h_c=6.3, z=6.3, θ=0, x=2.6, h_ho=2.6, material='cement', wall_load='max_normal_vert_wall', d_c=3.6,
                 d_ho=0.200, shape='circular')


plot_graphs.plot_graphs()
