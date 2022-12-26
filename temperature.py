from constants import Coefficients
from typing import Union

class Temperature(Coefficients):
    def __init__(self, in_temp=60, dt=1, DS_flow_rate=0.031) -> None:
        super().__init__(in_temp=in_temp, dt=dt, DS_flow_rate=DS_flow_rate)
        """
        Get important constants and topology
        """
        self.num_cells=50
        #initialize starting temperature at equilibrium, based of a gradient
        self.t_DS=[in_temp+(i*1.27) for i in range(self.num_cells)]
        self.t_DS_wall=[in_temp+(i*1.27) for i in range(self.num_cells)]
        self.t_AN=[in_temp+(i*1.27) for i in range(self.num_cells)]
        self.t_BW=[in_temp+(i*1.27) for i in range(self.num_cells)] 
        self.Ar=[self.t_DS,self.t_DS_wall,self.t_AN,self.t_BW]

        self.fahrenheit_to_kelvin(self.Ar)
        self.calc_temp()
            
    def fahrenheit_to_kelvin(self, temp: Union[int, float, list]):
        if type(temp) == list:
            if type(temp[0]) == list:
                for i in range(len(temp)):
                    temp[i] = [(((x-32)*5)/9)+273.15 for x in temp[i]]
                return temp
            elif type(temp[0]) == int or type(temp[0]) == float: 
                temp = [(((x-32)*5)/9)+273.15 for x in temp]
                return temp
        
        elif type(temp) == float or type(temp) == int:
            temp = (((temp-32)*5)/9)+273.15
            return temp

        else:
            raise ValueError("Fahrenheit to Kelvin function only accepts ints floats or lists")
    
    def kelvin_to_fahrenheit(self, temp: Union[int, float, list]):
        if type(temp) == list:
            if type(temp[0]) == list:
                for i in range(len(temp)):
                    temp[i] = [(((x-273.15)*9)/5)+32 for x in temp[i]]
                return temp
            elif type(temp[0]) == int or type(temp[0]) == float: 
                temp = [(((x-273.15)*9)/5)+32 for x in temp]
                return temp
        
        elif type(temp) == float or type(temp) == int:
            temp = (((temp-273.15)*9)/5)+32
            return temp

        else:
            raise ValueError("Kelvin to Fahrenheit function only accepts ints floats or lists")
    
    def calc_temp_DS(self, step):
        coef = self.coef_DS(step)

        if step == 0:
            t1 = self.Ar[0][0]
            t2 = self.fahrenheit_to_kelvin(60)
            t3 = self.Ar[1][0]

            t4 = ((t1*coef[0])+(t2*coef[1])+(t3*coef[2])+coef[3])*coef[4]

            self.Ar[0][step] = t4
            #print(self.kelvin_to_fahrenheit(t4),self.kelvin_to_fahrenheit(t1),self.kelvin_to_fahrenheit(t2),self.kelvin_to_fahrenheit(t3) )
        
        else:
            t1 = self.kelvin_to_fahrenheit(self.Ar[0][step])
            t2 = self.kelvin_to_fahrenheit(self.Ar[0][step-1])
            t3 = self.kelvin_to_fahrenheit(self.Ar[1][step])
            
            t4 = ((t1*coef[0])+(t2*coef[1])+(t3*coef[2])+coef[3])*coef[4]
            self.Ar[0][step] = self.fahrenheit_to_kelvin(t4)
           # if step == 1:
                #print(t4,coef,((t1*coef[0])+(t2*coef[1])+(t3*coef[2])+coef[3])* 0.02718483757768581,t1,t2,t3)
            
    def calc_temp_DS_wall(self, step):
        coef = self.coef_DS_wall(step)

        if step == 0:
            t1 = self.Ar[1][0]
            t2 = self.Ar[1][1]
            t3 = self.Ar[0][0]
            t4 = self.Ar[2][0]
            t5 = self.Ar[1][0]

            t6 = ((t1*coef[0])+(t2*coef[1])+(t3*coef[2])+(t4*coef[3])+(t5*coef[4]))*coef[5]
            self.Ar[1][0] = t6

        else: 
            t1 = self.Ar[1][step-1]
            if step < self.num_cells-1:
                t2 = self.Ar[1][step+1]
            else:
                t2 = self.Ar[1][step-1]
            t3 = self.Ar[0][step]
            t4 = self.Ar[2][step]
            t5 = self.Ar[1][step]

            t6 = ((t1*coef[0])+(t2*coef[1])+(t3*coef[2])+(t4*coef[3])+(t5*coef[4]))*coef[5]
            self.Ar[1][step] = t6

    def calc_temp_AN(self, step):
        coef = self.coef_AN(step)

        if step == 0:
            t1 = self.Ar[2][0]
            t2 = self.Ar[1][0]
            t3 = self.Ar[3][0]
            t4 = self.Ar[2][1]

            t5 = ((t1*coef[0])+(t2*coef[1])+(t3*coef[2])+(t4*coef[3])+coef[4])*coef[5]
            self.Ar[2][0] = t5
            #print(self.kelvin_to_fahrenheit(t5),coef,step)

        elif step == self.num_cells-1:
            self.Ar[2][self.num_cells-1] = self.Ar[0][self.num_cells-1]
            #print(self.Ar[2][step-1], self.Ar[1][step-1], step)

        else:
            t1 = self.Ar[2][step]
            if step < self.num_cells-1:
                t4 = self.Ar[2][step+1]
            else: 
                t4 = self.Ar[0][step]
            t2 = self.Ar[1][step]
            t3 = self.Ar[3][step]

            t5 = ((t1*coef[0])+(t2*coef[1])+(t3*coef[2])+(t4*coef[3])+coef[4])*coef[5]
            self.Ar[2][step] = t5
            

    def calc_temp(self):
        for i in range(self.num_cells):
            self.calc_temp_DS(i)
            self.calc_temp_DS_wall(i)
            self.calc_temp_AN(i)

