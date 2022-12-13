from topology.topology import Topology
import math

class Constants(Topology):
    """
    This class contains all the physical constants and
    physial properties of the well
    topology class should be an inout here 

    terminoly:
        -constants starting with (1,2,3) refer to (DS and AN fluid, DS wall, and bore wall)
    """
    def __init__(self,in_temp=60,dt=1,DS_flow_rate=0.05257) -> None:
        """
        Calculate and get physical values such as:
        -radii for all cylinders
        -dz
        -max depth of well
        """

        #parsing topology data 
        bore_well=Topology()
        bore_well.populate_radii()
        bore_well.calc_cross_sections()
        self.topology=bore_well.topology
        self.dt=dt
        self.dz=bore_well.dz
        self.height=bore_well.height
        self.num_cells=bore_well.num_cells

        #inputs that should come from pressure or density calcs
        self.den_DS = 805
        self.den_DS_wall = 7000
        self.den_An = 805

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

        #possible future variables
        self.q = DS_flow_rate
        self.GG = 0.0127
        self.supTemp = in_temp*255.928
        self.Q1 = 1
        self.Q2 = 1

class Coefficients(Constants):
    """
    Contains all the coefficients for temperature
    soon it will include some 
    """
    def __init__(self) -> None:
        self.constants = Constants()

    def coef_DS(self,cell):
        inner_radius=self.get_inner_radius(cell)

        A = (self.constants.den_DS*self.constants.c1)/self.constants.dt
        B = (self.constants.den_DS*self.constants.c1*self.constants.q)*(math.pi*(inner_radius**2)/self.constants.dz)
        C = (2*self.constants.h1)*inner_radius
        D = self.constants.Q1*(2*math.pi*(inner_radius**2))
        E = 1/(A+B+C) 
        return [A,B,C,D,E]

    def coef_DS_wall(self,cell): #someting is wrong here
        inner_radius=self.get_inner_radius(cell)
        outer_radius=self.get_outer_radius(cell)

        A1 = (self.constants.l2)/(self.constants.dz*(self.constants.dz/2))
        B1 = (self.constants.l2)/(self.constants.dz*(self.constants.dz/2))
        C1 = (2*inner_radius*self.constants.h1)*((outer_radius**2)-(inner_radius**2))
        E1 = (2*outer_radius*self.constants.h2)*((outer_radius**2)-(inner_radius**2))
        F1 = (self.constants.den_DS_wall*self.constants.c2)/self.constants.dt
        G1 = 1/(A1+B1+C1+E1+F1)
        return [A1,B1,C1,E1,F1,G1]

    def coef_AN(self,cell):
        outer_radius=self.get_outer_radius(cell)
        bore_radius=self.get_bore_radius(cell)

        A = (self.constants.c3*self.constants.den_An)/self.constants.dt
        B = (2*self.constants.h2)*outer_radius
        C = (2*self.constants.h3)*bore_radius
        D = (self.constants.c3*self.constants.den_An*self.constants.q)*(math.pi*((bore_radius**2)-(outer_radius**2)))/self.constants.dz
        E = self.constants.Q2*(math.pi*((bore_radius**2)-(outer_radius**2)))
        F = 1/(A+B+C+D)
        return [A,B,C,D,E,F]

    
    """
    Access topology data when needed
    """
    def get_inner_radius(self,cell):
       # print(cell,self.constants.topology)
        return self.constants.topology['inner radii'][cell]
    
    def get_outer_radius(self,cell):
        return self.constants.topology['outer radii'][cell]

    def get_bore_radius(self,cell):
        return self.constants.topology['bore radii'][cell]

    def get_inner_area(self,cell):
        return self.constants.topology['DS areas'][cell]

    def get_inner_area(self,cell):
        return self.constants.topology['AN areas'][cell]
    
    def get_inner_area(self,cell):
        return self.constants.topology['pipe areas'][cell]

