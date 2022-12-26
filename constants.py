from topology import Topology
import math

class Constants(Topology):
    """
    This class contains all the physical constants and
    physial properties of the well
    topology class should be an inout here 

    terminoly:
        -constants starting with (1,2,3) refer to (DS and AN fluid, DS wall, and bore wall)
    """
    def __init__(self, inner_radii_boundaries=[(0,4),(3500,3.5),(4500,2.5)], outer_radii_boundaries=[(0,5),(3500,5.5),(4500,6)], bore_wall_boundaries=[(0,10),(3500,8.5)], height=5000.0, dz=100, pipe_roughness=0.0015,in_temp=60,dt=1,DS_flow_rate=.0315451):
        super().__init__(inner_radii_boundaries, outer_radii_boundaries, bore_wall_boundaries, height, dz, pipe_roughness)
        """
        Calculate and get physical values such as:
        -radii for all cylinders
        -dz
        -max depth of well
        """
        self.dt=dt

        #inputs that should come from pressure or density calcs
        self.den_DS = 805    #[kg/m^3]
        self.den_DS_wall = 7000
        self.den_AN = 805

        #physical constants, classified by:
            #convection factor
            #thermal conductivity
            #specific heat capacity 
        self.h1 = 3
        self.h2 = 13
        self.h3 = 3
        self.hef = 0.4

        self.l1 = 8.1
        self.l2 = 8.09
        self.l3 = 8.2

        self.c1 = 16
        self.c2 = 30 #44
        self.c3 = 20

        self.grav=9.81
        #self.gas_constant=8.31145

        #possible future variables
        self.q = DS_flow_rate
        self.GG = 0.0127
        self.supTemp = in_temp*255.928
        self.Q1 = 1
        self.Q2 = 1
        self.mat_viscosity=0.015 #[kg/m*s]
        self.surface_tension=0.034

class Coefficients(Constants):
    """
    Contains all the coefficients for temperature
    soon it will include some 
    """
    def __init__(self,in_temp=60,dt=1,DS_flow_rate=0.05257) -> None:
        super().__init__(in_temp=in_temp, dt=dt, DS_flow_rate=DS_flow_rate)

    def coef_DS(self,cell):
        inner_radius = self.inner_radii[cell]

        A = (self.den_DS*self.c1)/self.dt
        B = (self.den_DS*self.c1*self.q)*(math.pi*(inner_radius**2)/self.dz)
        C = (2*self.h1)*inner_radius
        D = self.Q1*(2*math.pi*(inner_radius**2))
        E = 1/(A+B+C) 
        return [A,B,C,D,E]

    def coef_DS_wall(self,cell): #someting is wrong here
        inner_radius=self.inner_radii[cell]
        outer_radius=self.outer_radii[cell]

        A1 = (self.l2)/(self.dz*(self.dz/2))
        B1 = (self.l2)/(self.dz*(self.dz/2))
        C1 = (2*inner_radius*self.h1)*((outer_radius**2)-(inner_radius**2))
        E1 = (2*outer_radius*self.h2)*((outer_radius**2)-(inner_radius**2))
        F1 = (self.den_DS_wall*self.c2)/self.dt
        G1 = 1/(A1+B1+C1+E1+F1)
        return [A1,B1,C1,E1,F1,G1]

    def coef_AN(self,cell):
        outer_radius=self.outer_radii[cell]
        bore_radius=self.bore_radii[cell]

        A = (self.c3*self.den_AN)/self.dt
        B = (2*self.h2)*outer_radius
        C = (2*self.h3)*bore_radius
        D = (self.c3*self.den_AN*self.q)*(math.pi*((bore_radius**2)-(outer_radius**2)))/self.dz
        E = self.Q2*(math.pi*((bore_radius**2)-(outer_radius**2)))
        F = 1/(A+B+C+D)
        return [A,B,C,D,E,F]

