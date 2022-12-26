from constants import Constants
import matplotlib.pyplot as plt

class Velocity(Constants): 
    """
    We are taking the inlet pressure of the pump 
    and taking that to be the velocity
    this is equivalent to q in Constants 
    """
    def __init__(self, in_temp=60, dt=1, DS_flow_rate=0.05257, DS_gas_flow_rate=0.031) -> None:
        super().__init__(in_temp=in_temp, dt=dt, DS_flow_rate=DS_flow_rate)
        self.q_gas=DS_gas_flow_rate
        self.liquid_sup_velocities={'DS superficial velocities': [], 'AN superficial velocities': []}
        self.gas_sup_velocities={'DS superficial velocities': [], 'AN superficial velocities': []}

        self.calc_sup_velocities()
        #print(self.liquid_sup_velocities)

    def calc_DS_liquid_sup_velocity(self,cell):              #does not take into account friction
        area=self.topology['DS areas'][cell]*0.00064516
        if area <= 0:
            raise Exception('Area cannot be zero or negative')
        velocity=self.q/area
        self.liquid_sup_velocities['DS superficial velocities'].append(velocity)
    
    def calc_AN_liquid_sup_velocity(self,cell):              #does not take into account friction or gravity or 
        area=self.topology['AN areas'][cell]*0.00064516
        if area <= 0:
            raise Exception('Area cannot be zero or negative')
        velocity=self.q/area
        self.liquid_sup_velocities['AN superficial velocities'].append(velocity)

    def calc_DS_gas_sup_velocity(self,cell):              #does not take into account friction
        area=self.topology['DS areas'][cell]*0.00064516
        if area <= 0:
            raise Exception('Area cannot be zero or negative')
        velocity=self.q_gas/area
        self.gas_sup_velocities['DS superficial velocities'].append(velocity)
    
    def calc_AN_gas_sup_velocity(self,cell):              #does not take into account friction or gravity or 
        area=self.topology['AN areas'][cell]*0.00064516
        if area <= 0:
            raise Exception('Area cannot be zero or negative')
        velocity=self.q_gas/area
        self.gas_sup_velocities['AN superficial velocities'].append(velocity)
    
    def calc_sup_velocities(self):
        for i in range(self.num_cells):
            self.calc_DS_liquid_sup_velocity(i)
            self.calc_DS_gas_sup_velocity(i)
            self.calc_AN_liquid_sup_velocity(i)
            self.calc_AN_gas_sup_velocity(i)    
