import math
from math import pi as π
from dataclasses import dataclass, field

@dataclass
class SiloCalcs:
    ''' Assumptions:
    1. Silos are filled and emptied from centre
    2. Silos are flat-topped
    '''

    # input - bulk solid  parameters
    material: str

    # input - silo dimensions 
    d_c: float  # silo inside critical dimension (m)
    d_ho: float  # silo hopper outlet diameter (m)
    h_ho: float  # hopper height transition to outlet (m)
    shape: str  # shape of silo - 'circular' or 'square'
    h_c: float  # height of vertical walled segment (m)

    # input - wall loading assessment
    wall_load: str
    θ: float = 0  # patch load circumferential coordinate (rad)
    z: float = -1  # depth of pressure analyse below equivalent surface
    x: float = -1  # height above hopper apex

    # calculated values from inputs
    def __post_init__(self):
        if self.material == 'cement':
            self.γ_l = 13.0
            self.γ_u = 16.0
            self.ψ_r = 36
            self.ψ_im = 30
            self.a_ψ = 1.22
            self.K_m = 0.54
            self.a_K = 1.20
            self.μ_m = 0.46
            self.a_μ = 1.07
            self.C_op = 0.5

        elif self.material == 'flyash':
            self.γ_l = 8.0
            self.γ_u = 15.0
            self.ψ_r = 41
            self.ψ_im = 35
            self.a_ψ = 1.16
            self.K_m = 0.46
            self.a_K = 1.20
            self.μ_m = 0.62
            self.a_μ = 1.07
            self.C_op = 0.5

        elif self.material == 'sand':
            self.γ_l = 14.0
            self.γ_u = 16.0
            self.ψ_r = 39
            self.ψ_im = 36
            self.a_ψ = 1.09
            self.K_m = 0.45
            self.a_K = 1.11
            self.μ_m = 0.48
            self.a_μ = 1.16
            self.C_op = 0.4
        else:
            print("Error: Material not defined")

        # characteristic value to be adopted - Table 3.1
        self.γ = self.γ_u
        if self.wall_load == 'max_normal_vert_wall':
            self.μ = self.μ_m / self.a_μ
            self.K = self.K_m * self.a_K
            self.ψ = self.ψ_im / self.a_ψ
        elif self.wall_load == 'max_friction_vert_wall':
            self.μ = self.μ_m * self.a_μ
            self.K = self.K_m * self.a_K
            self.ψ = self.ψ_im / self.a_ψ
        elif self.wall_load == 'max_vert_load_hopper':
            self.μ = self.μ_m / self.a_μ
            self.K = self.K_m / self.a_K
            self.ψ = self.ψ_im * self.a_ψ
        elif self.wall_load == 'max_hopper_filling':
            self.μ = self.μ_m / self.a_μ
            self.K = self.K_m / self.a_K
            self.ψ = self.ψ_im / self.a_ψ
        elif self.wall_load == 'max_hopper_discharge':
            self.μ = self.μ_m / self.a_μ
            self.K = self.K_m * self.a_K
            self.ψ = self.ψ_im * self.a_ψ
        else:
            print("Error: wall loading assessment not defined")

        self.r = self.d_c / 2  # equivalent radius of silo (m)
        self.h_h = self.h_ho + (self.h_ho * self.d_ho) / (self.d_c - self.d_ho)  # hopper apex height (m) - similar triangles
        self.β = math.atan(self.r / self.h_h)  # hopper apex half angle (rad)

        self.e_f = 0.25 * self.d_c  # maximum eccentricity of surface pile during filling
        self.s = 0.2 * self.d_c  # height of patch load zone

        #convert to radians
        self.ψ_r = math.radians(self.ψ_r)
        self.ψ = math.radians(self.ψ)

        if self.shape == 'circular':
            self.A = π * self.r**2
            self.U = π * self.d_c
        elif self.shape == 'square':
            self.A = self.d_c**2
            self.U = self.d_c**2
        else:
            print("Error: Silo shape not defined")

        if self.z < 0:
            self.z = self.h_c

        if self.x < 0:
            self.x = self.h_ho

    def plot_values(self):
        return [self.h_c, self.h_ho]

    def usable_volume(self):
        if self.shape == 'circular':
            # Hopper geometry calculations
            h_h = self.r / math.tan(self.β)  # hopper height
            V_h = (1/3) * π * h_h * self.r**2  # hopper volume
            # Bulk solids top pile geometry calculations
            h_tp = self.r * math.tan(self.ψ_r)  # height of pile on top of silo due to internal friction
            V_tp = (1/3) * π * h_tp * self.r**2  # volume of pile
            # Silo usable volume
            h_sf = 0.5  # silo height of saftey gap below ceiling (m)
            V_vw = self.A * (self.h_c - h_sf - h_tp)  # volume of vertical wall section
            V_t = V_h + V_tp + V_vw  # total useable volume of silo
            return V_t

        elif self.shape == 'square':
            A = self.d_c**2
            # Hopper geometry calculations
            h_h = self.r / math.tan(self.β)  # hopper height
            V_h = (1/3) * h_h * self.d_c**2 
            # Bulk solids top pile geometry calculations
            h_tp = self.r * math.tan(self.ψ_r)  # height of pile on top of silo due to internal friction
            V_tp = (1/3) * h_tp * self.d_c**2  # volume of pile
            # Silo usable volume
            h_sf = 0.5  # silo height of saftey gap below ceiling (m)
            V_vw = self.A * (self.h_c - h_sf - h_tp)  # volume of vertical wall section
            V_t = V_h + V_tp + V_vw  # total usable volume of silo
            return V_t

    def silo_capacity(self):
        V_t = self.usable_volume()
        min_capacity = V_t * self.γ_l / 9.81
        max_load = V_t * self.γ_u / 9.81
        d = {
            'min_capacity': min_capacity,
            'max_load': max_load
        }
        return d

    def f_z0(self):  # Janssen characteristic depth
        z0 = (1 / (self.K * self.μ)) * (self.A / self.U)  # Eq.5.5
        return z0

    def f_Y_Jz(self):  # Janssen pressure depth variation
        Y_Jz = 1 - math.exp(- self.z / self.f_z0())  # Eq.5.6
        return Y_Jz

    def f_p_h0(self):  # asymptotic horizontal pressure at great depth
        p_h0 = self.γ * self.K * self.f_z0()  # Eq.5.4
        return p_h0

    def f_E(self):
        E = 2 * (self.e_f / self.d_c)  # Eq.5.10
        return E

    def WallFillingLoad(self):
        p_h0 = self.f_p_h0()
        Y_Jz = self.f_Y_Jz()
        p_hf = p_h0 * Y_Jz  # horizontal pressure Eq.5.1
        p_wf = self.μ * p_h0 * Y_Jz  # wall friction traction Eq.5.2
        p_vf = (1 / self.K) * p_h0 * Y_Jz  # vertical pressure Eq.5.3
        n_zSk = self.μ * p_h0 * (self.z - self.f_z0() * self.f_Y_Jz())  # vertical force in wall per unit of perimeter
        d = {
            'p_hf': p_hf,
            'p_wf': p_wf,
            'p_vf': p_vf,
            'n_zSk': n_zSk
        }
        return d

    def FillingPatchLoad(self):
        z_p = min(self.f_z0(), 0.5 * self.h_c)  # depth of patch zone below eq. surface Eq.5.16
        E = self.f_E()
        C_pf = 0.21*self.C_op*(1+2*E**2)*(1-math.exp(-1.5*((self.h_c/self.d_c) - 1)))
        if C_pf < 0:
            C_pf = 0
        p_pf = C_pf * self.WallFillingLoad()['p_hf']  # filling patch pressure magnitude Eq.5.8
        p_pfs = p_pf * math.cos(self.θ)  # circumferential patch pressure variation Eq.5.14
        F_pf = (π / 2) * self.s * self.d_c * p_pf  # total horizontal force due to filling patch pressure
        if self.shape == 'circular':
            d = {
                'p_pf': p_pf,
                'p_pfs': p_pfs,
                'F_pf': F_pf,  # total horizontal load
                'z_p': z_p  # depth of patch load below equivalent surface
            }
        elif self.shape == 'sqaure':
            d = {
                'p_pf_nc': 0.36 * p_pf  # non circular filling patch band load Eq.5.17
            }
        return d

    def WallDischargeLoad(self):
        C_h = 1.15  # horizontal discharge factor Eq.5.21
        C_w = 1.10  # wall friction discharge factor Eq.5.22
        p_he = C_h * self.WallFillingLoad()['p_hf']  # Eq.5.18
        p_we = C_w * self.WallFillingLoad()['p_wf']  # Eq.5.19
        n_zSk = C_w * self.WallFillingLoad()['n_zSk']
        d = {
            'p_he': p_he,
            'p_we': p_we,
            'n_zSk': n_zSk
        }
        return d

    def DischargePatchLoad(self):
        z_p = min(self.f_z0(), 0.5 * self.h_c)  # depth of patch zone below eq. surface Eq.5.16
        E = self.f_E()
        C_pe = 0.42*self.C_op*(1+2*E**2)*(1-math.exp(-1.5*((self.h_c/self.d_c) - 1)))  # Eq.5.28
        p_pe = C_pe * self.WallDischargeLoad()['p_he']  # p_he is local value at height that patch load is applied Eq.5.27
        p_pes = p_pe * math.cos(self.θ)
        F_pe = (π / 2) * self.s * self.d_c * p_pe  # total horizontal force due to filling patch pressure
        if self.shape == 'circular':
            d = {
                'p_pe': p_pe,
                'p_pes': p_pes,
                'F_pe': F_pe,
                'z_p': z_p
            }
        elif self.shape == 'square':
            d = {
                'p_pe_nc': 0.36 * p_pe  # non circular filling patch band load Eq.5.37
            }
        return d

    def HopperFillingLoad(self):
        if (math.tan(self.β) > (1-self.K)/(2 * self.μ)):
            self.μ = (1-self.K) / (2 * math.tan(self.β))  # shallow hopper Eq.6.26
        if self.z != self.h_c:
            print("z is not equal to the transition and results are not valid!")
        C_b = 1.0  # bottom load magnifier Eq.6.3
        p_vft = C_b * self.WallFillingLoad()['p_vf']
        b = 0.2  # empirical coefficient §6.3.2.
        F_f = 1 - b / (1 + math.tan(self.β) / (self.μ))
        S = 2  # value for conical hoppers Eq.6.9
        n = S * (1 - b) * self.μ / math.tan(self.β)  # n = S*(F_f * self.μ / math.tan(self.β) + F_f) - 2
        p_v = ((self.γ * self.h_h )/(n-1))*((self.x/self.h_h)-(self.x/self.h_h)**n)+p_vft*(self.x/self.h_h)**n  # mean vertical stress in solid at height x above apex
        p_nf = F_f * p_v  # normal pressure in hopper during filling and full load
        p_tf = self.μ * F_f * p_v  # frictional traction in hopper at x along hopper wall
        d = {
            'p_v': p_v,
            'p_nf': p_nf,
            'p_tf': p_tf
        }
        return d

    def HopperDischargeLoad(self):
        if (math.tan(self.β) > (1-self.K)/(2 * self.μ)):
            self.μ = (1-self.K) / (2 * math.tan(self.β))  # shallow hopper Eq.6.26
        if self.z != self.h_c:
            print("z is not equal to the transition and results are not valid!")
        C_b = 1.0  # bottom load magnifier Eq.6.3
        p_vft = C_b * self.WallFillingLoad()['p_vf']
        ϕ_wh = math.atan(self.μ)
        ε = ϕ_wh + math.asin(math.sin(ϕ_wh) / math.sin(self.ψ))
        F_e = (1 + math.sin(self.ψ) * math.cos(ε)) / (1 - math.sin(self.ψ) * math.cos(2 * self.β + ε))
        S = 2  # value for conical hoppers Eq.6.9
        n = S*(F_e * self.μ / math.tan(self.β) + F_e) - 2
        p_v = ((self.γ * self.h_h )/(n-1))*((self.x/self.h_h)-(self.x/self.h_h)**n)+p_vft*(self.x/self.h_h)**n  # mean vertical stress in solid at height x above apex
        p_ne = F_e * p_v  # normal pressure in hopper during filling and full load
        p_te = self.μ * F_e * p_v  # frictional traction in hopper at x along hopper wall
        d = {
            'p_v': p_v,
            'p_ne': p_ne,
            'p_te': p_te
        }
        return d

if __name__ == "__main__":
    print("test")