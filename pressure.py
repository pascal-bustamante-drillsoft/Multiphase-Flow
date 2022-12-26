from velocity import Velocity
#from flow import Flow
import math

class Pressure(Velocity):
    def __init__(self, in_temp=60, dt=1, DS_flow_rate=0.031) -> None:   #check dimensions and write comments
        super().__init__(in_temp=in_temp, dt=dt, DS_flow_rate=DS_flow_rate)
        self.static_pressure=[]
        self.calc_static_pressure()
        #print(self.static_pressure)


    def calc_local_static_pressure(self,cell):
        pressure=self.dz*0.3048*self.den_DS*9.81
        if cell == 0:
            if len(self.static_pressure)==0:
                self.static_pressure.append(pressure)
            else:
                raise Exception('Indexing issue in static pressure calculation')
        else:
            if cell == len(self.static_pressure):
                self.static_pressure.append(self.static_pressure[cell-1]+pressure)
            else:
                raise Exception('Indexing issue in static pressure calculation')
    
    def calc_static_pressure(self):
        for i in range(self.num_cells):
            self.calc_local_static_pressure(i)

