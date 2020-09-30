import numpy as np
import matplotlib.pyplot as plt
from silo_calcs import SiloCalcs

# plot variables
plot_material = 'cement'

plot_h_c = 6.3  # height of vertical wall
plot_h_ho = 2.6  # height of hopper to outlet

plot_d_c = 3.6  # critical dimension of silo
plot_d_ho = 0.2  # size of hopper outlet

plot_shape = 'circular'


def plot_graphs():

    x_i_lst = np.linspace(0, plot_h_ho, 100)  # x coordinates for plotting graph
    z_i_lst = np.linspace(0, plot_h_c, 100)  # z coordinates for plotting graph

    # wall filling loads
    wall_max_pn_f = []
    wall_max_pft_f = []
    wall_max_pv_f = []

    # wall discharge loads
    wall_max_pn_e = []
    wall_max_pft_e = []
    wall_max_pv_e = []

    # hopper filling loads
    hopper_max_pv_f = []
    hopper_max_pn_f = []
    hopper_max_pt_f = []

    # hopper discharge loads
    hopper_max_pv_e = []
    hopper_max_pn_e = []
    hopper_max_pt_e = []



    for idz, z_i in enumerate(z_i_lst):
        # normal force on vertical wall
        silo_p_hf = SiloCalcs(z=z_i, material=plot_material, d_c=plot_d_c, d_ho=plot_d_ho, h_ho=plot_h_ho,
                              shape=plot_shape, h_c=plot_h_c, wall_load='max_normal_vert_wall')
        wall_max_pn_f.append(silo_p_hf.WallFillingLoad()['p_hf'])
        wall_max_pn_e.append(silo_p_hf.WallDischargeLoad()['p_he'])
        # friction traction force on vertical wall
        silo_p_wf = SiloCalcs(z=z_i, material=plot_material, d_c=plot_d_c, d_ho=plot_d_ho, h_ho=plot_h_ho,
                              shape=plot_shape, h_c=plot_h_c, wall_load='max_friction_vert_wall')
        wall_max_pft_f.append(silo_p_wf.WallFillingLoad()['p_wf'])
        wall_max_pft_e.append(silo_p_wf.WallDischargeLoad()['p_we'])
        # vertical force on hopper
        silo_p_vf = SiloCalcs(z=z_i, material=plot_material, d_c=plot_d_c, d_ho=plot_d_ho, h_ho=plot_h_ho,
                              shape=plot_shape, h_c=plot_h_c, wall_load='max_vert_load_hopper')
        wall_max_pv_f.append(silo_p_vf.WallFillingLoad()['p_vf'])
        wall_max_pv_e.append(silo_p_vf.WallFillingLoad()['p_vf'])

    for idx, x_i in enumerate(x_i_lst):
        hopper_filling = SiloCalcs(x=x_i, material=plot_material, d_c=plot_d_c, d_ho=plot_d_ho, h_ho=plot_h_ho,
                              shape=plot_shape, h_c=plot_h_c, wall_load='max_hopper_filling')
        hopper_max_pv_f.append(hopper_filling.HopperFillingLoad()['p_v'])
        hopper_max_pn_f.append(hopper_filling.HopperFillingLoad()['p_nf'])
        hopper_max_pt_f.append(hopper_filling.HopperFillingLoad()['p_tf'])

        hopper_discharge = SiloCalcs(x=x_i, material=plot_material, d_c=plot_d_c, d_ho=plot_d_ho, h_ho=plot_h_ho,
                              shape=plot_shape, h_c=plot_h_c, wall_load='max_hopper_discharge')
        hopper_max_pv_e.append(hopper_discharge.HopperDischargeLoad()['p_v'])
        hopper_max_pn_e.append(hopper_discharge.HopperDischargeLoad()['p_ne'])
        hopper_max_pt_e.append(hopper_discharge.HopperDischargeLoad()['p_te'])

    # rearrange lists to use common coordinates
    z_p_lst_i = z_i_lst + plot_h_ho
    z_p_lst = np.concatenate((x_i_lst, z_p_lst_i))

    # rearrange filling pressures for plot1
    plot_max_pn_f = hopper_max_pn_f
    plot_max_pn_f.extend(wall_max_pn_f[::-1])

    plot_max_pft_f = hopper_max_pt_f
    plot_max_pft_f.extend(wall_max_pft_f[::-1])

    plot_max_pv_f = hopper_max_pv_f
    plot_max_pv_f.extend(wall_max_pv_f[::-1])

    plot1 = plt.figure(1) 
    plt.plot(plot_max_pn_f, z_p_lst, label="p_nf")
    plt.plot(plot_max_pft_f, z_p_lst, label="p_tf")
    plt.plot(plot_max_pv_f, z_p_lst, label="p_vf")

    plt.xlabel('Pressure (kN/m3) (kPa)')
    plt.ylabel('Height above hopper apex (m)')
    plt.title('Symmetrical filling load')
    plt.legend()

    # rearrange discharge pressures for plot2
    plot_max_pn_e = hopper_max_pn_e
    plot_max_pn_e.extend(wall_max_pn_e[::-1])

    plot_max_pft_e = hopper_max_pt_e
    plot_max_pft_e.extend(wall_max_pft_e[::-1])

    plot_max_pv_e = hopper_max_pv_e
    plot_max_pv_e.extend(wall_max_pv_e[::-1])

    plot2 = plt.figure(2) 
    plt.plot(plot_max_pn_e, z_p_lst, label="p_ne")
    plt.plot(plot_max_pft_e, z_p_lst, label="p_te")
    plt.plot(plot_max_pv_e, z_p_lst, label="p_ve")

    plt.xlabel('Pressure (kN/m3) (kPa)')
    plt.ylabel('Height above hopper apex (m)')
    plt.title('Symmetrical discharge load')
    plt.legend()

    plt.show()
