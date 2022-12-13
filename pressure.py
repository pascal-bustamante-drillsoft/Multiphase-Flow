from grid import Topology
from constants import Constants


class Pressure(Constants):
    def __init__(self):
        self.topology=Topology()
        self.grid=self.topology.grid
        self.constants=Constants()
        self.pressure_DS=[] 
        self.pressure_AN=[]
        self.cross_sectional_area_DS = self.calc_cross_sectional_area_DS()
        
        inner_radii=self.grid['inner_radii']

    def calc_cross_sectional_area_DS(self):