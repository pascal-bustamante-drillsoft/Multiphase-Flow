from temperature import Temperature
from pressure import Pressure 
import matplotlib.pyplot as plt
import math

class Gas(Temperature,Pressure):
    def __init__(self, in_temp=60, dt=1, DS_flow_rate=0.031) -> None:
        super().__init__(in_temp=in_temp, dt=dt, DS_flow_rate=DS_flow_rate)
        """
        Specific gravity of natural gas based on the ratio of the molecular weight of natural gas and air
        """
        self.gas_constant=8.31145
        self.specific_gravity=0.5
        self.z_factors={'DS z factors':[],'AN z factors':[]}
        self.volumes={'DS volumes':[],'AN volumes':[]}
        self.gas_densities={'DS densities':[],'AN densities':[]}
        self.densities=[]
        self.viscosities={'DS viscosities':[],'AN viscosities':[]}

        self.calc_z_factors()
        self.calc_volumes()
        self.calc_densities()
        self.calc_viscosities()
        #print(self.liquid_sup_velocities)
    
    def calc_local_z_factor(self,cell):
        """
        Based of off Coker, A K (1993) Program calculates Z-factor for natural gas. (for natural gas only!)
        pressure cannot exceede 5000psia for his work to apply
        """
        #local_pressure=self.static_pressure[cell]
        ##TEMP!!
        pressure=self.static_pressure
        #pressure=[240702.08400000003, 481404.16800000006, 722106.2520000001, 962808.3360000001, 1203510.4200000002, 1444212.5040000002, 1684914.5880000002, 1925616.6720000003, 2166318.756, 2407020.84, 2647722.9239999996, 2888425.0079999994, 3129127.0919999992, 3369829.175999999, 3610531.259999999, 3851233.3439999986, 4091935.4279999984, 4332637.511999998, 4573339.595999998, 4814041.679999998, 5054743.763999998, 5295445.847999997, 5536147.931999997, 5776850.015999997, 6017552.099999997, 6258254.183999997, 6498956.267999996, 6739658.351999996, 6980360.435999996, 7221062.519999996, 7461764.603999996, 7702466.687999995, 7943168.771999995, 8183870.855999995, 8424572.939999996, 8665275.023999996, 8905977.107999997, 9146679.191999998, 9387381.275999999, 9628083.36, 9868785.444, 10109487.528, 10350189.612000002, 10590891.696000002, 10831593.780000003, 11072295.864000004, 11312997.948000005, 11553700.032000005, 11794402.116000006, 12035104.200000007]
        local_pressure=(pressure[cell]*0.000145038)/1000


        local_DS_temp=self.t_DS[cell]
        local_AN_temp=self.t_AN[cell]

        """
        Calc needed coefs
        """
        F1=(local_pressure)*(0.251*self.specific_gravity-0.15)-0.202*self.specific_gravity+1.106
        F2_DS=1.4*math.exp(-0.0054*(local_DS_temp))
        F3=0.001946*(local_pressure**5)-0.027635*(local_pressure**4)+0.136315*(local_pressure**3)-0.238489*(local_pressure**2)+0.105168*(local_pressure)
        F4=(0.154-0.152*self.specific_gravity)*(local_pressure**(3.18*self.specific_gravity-1))*math.exp(-0.5*local_pressure)-0.02
        F5=0.35*((0.6-self.specific_gravity)*math.exp(-1.039*((local_pressure-1.8)**2)))
        F6_DS=1+((344000000*local_pressure*(10**(1.785*self.specific_gravity)))/((local_DS_temp+460)**3.825))

        F2_AN=1.4*math.exp(-0.0054*(local_AN_temp))
        F6_AN=1+((344000000*local_pressure*(10**(1.785*self.specific_gravity)))/((local_AN_temp+460)**3.825))

#        F6_AN=1+((344000000*local_pressure*(math.pow(10,1.785*self.specific_gravity)))/(math.pow(local_AN_temp,3.825)))

        local_z_factor_DS=F1*((1/F6_DS)+(F2_DS*F3))+F4+F5
        local_z_factor_AN=F1*((1/F6_AN)+(F2_AN*F3))+F4+F5
       # print(local_z_factor_DS,F1,F2_DS,F3,F4,F5,F6_DS,local_DS_temp,local_pressure)

        self.z_factors['DS z factors'].append(local_z_factor_DS)
        self.z_factors['AN z factors'].append(local_z_factor_AN)

    def calc_local_volume_per_mole(self,cell):
        '''
        Using the ideal gas law and Z-compresability factor
        we calculate the volume a mole of natural gas occupies
        the units are as following:
        -[P] = [atm]
        -[V] = [L]
        -[T] = [K]
        -[R] = [J/mol*K]
        -[n] = [mol]
        -Z is dimensionless
        '''
        pressure=self.static_pressure
        #pressure=[240702.08400000003, 481404.16800000006, 722106.2520000001, 962808.3360000001, 1203510.4200000002, 1444212.5040000002, 1684914.5880000002, 1925616.6720000003, 2166318.756, 2407020.84, 2647722.9239999996, 2888425.0079999994, 3129127.0919999992, 3369829.175999999, 3610531.259999999, 3851233.3439999986, 4091935.4279999984, 4332637.511999998, 4573339.595999998, 4814041.679999998, 5054743.763999998, 5295445.847999997, 5536147.931999997, 5776850.015999997, 6017552.099999997, 6258254.183999997, 6498956.267999996, 6739658.351999996, 6980360.435999996, 7221062.519999996, 7461764.603999996, 7702466.687999995, 7943168.771999995, 8183870.855999995, 8424572.939999996, 8665275.023999996, 8905977.107999997, 9146679.191999998, 9387381.275999999, 9628083.36, 9868785.444, 10109487.528, 10350189.612000002, 10590891.696000002, 10831593.780000003, 11072295.864000004, 11312997.948000005, 11553700.032000005, 11794402.116000006, 12035104.200000007]
        local_pressure=pressure[cell]*9.86923e-6

        local_DS_temp=self.fahrenheit_to_kelvin(self.t_DS[cell])
        local_AN_temp=self.fahrenheit_to_kelvin(self.t_AN[cell])

        local_z_factor_DS=self.z_factors['DS z factors'][cell]
        local_z_factor_AN=self.z_factors['AN z factors'][cell]

        local_volume_DS=(self.gas_constant*local_DS_temp*local_DS_temp*local_z_factor_DS)/local_pressure
        local_volume_AN=(self.gas_constant*local_AN_temp*local_AN_temp*local_z_factor_AN)/local_pressure

        self.volumes['DS volumes'].append(local_volume_DS*0.001)
        self.volumes['AN volumes'].append(local_volume_AN*0.001)

    def calc_local_gas_density(self,cell):
        pressure=self.static_pressure
        #pressure=[240702.08400000003, 481404.16800000006, 722106.2520000001, 962808.3360000001, 1203510.4200000002, 1444212.5040000002, 1684914.5880000002, 1925616.6720000003, 2166318.756, 2407020.84, 2647722.9239999996, 2888425.0079999994, 3129127.0919999992, 3369829.175999999, 3610531.259999999, 3851233.3439999986, 4091935.4279999984, 4332637.511999998, 4573339.595999998, 4814041.679999998, 5054743.763999998, 5295445.847999997, 5536147.931999997, 5776850.015999997, 6017552.099999997, 6258254.183999997, 6498956.267999996, 6739658.351999996, 6980360.435999996, 7221062.519999996, 7461764.603999996, 7702466.687999995, 7943168.771999995, 8183870.855999995, 8424572.939999996, 8665275.023999996, 8905977.107999997, 9146679.191999998, 9387381.275999999, 9628083.36, 9868785.444, 10109487.528, 10350189.612000002, 10590891.696000002, 10831593.780000003, 11072295.864000004, 11312997.948000005, 11553700.032000005, 11794402.116000006, 12035104.200000007]
        local_pressure=pressure[cell]*1e-6 #Mpa

        local_DS_temp=self.fahrenheit_to_kelvin(self.t_DS[cell])
        local_AN_temp=self.fahrenheit_to_kelvin(self.t_AN[cell])

        local_z_factor_DS=self.z_factors['DS z factors'][cell]
        local_z_factor_AN=self.z_factors['AN z factors'][cell]

        local_gas_density_DS=(3484.4*local_pressure*self.specific_gravity)/(local_z_factor_DS*local_DS_temp)
        local_gas_density_AN=(3484.4*local_pressure*self.specific_gravity)/(local_z_factor_AN*local_AN_temp)

        self.gas_densities['DS densities'].append(local_gas_density_DS)
        self.gas_densities['AN densities'].append(local_gas_density_AN)

    def calc_local_gas_viscosity(self,cell):  
        """
        Using Lee et al. (1966) correlation
        """
        local_DS_temp=self.t_DS[cell]+ 460
        local_AN_temp=self.t_AN[cell]+ 460
        local_DS_density=(self.gas_densities['DS densities'][cell])/1000
        local_AN_density=(self.gas_densities['AN densities'][cell])/1000
        avg_molecular_weight=19          #look up a good value and move them into an input based of the materials class
        """
        Calc coefs
        """
        #K_DS=(((7.77+0.0063*avg_molecular_weight)*(local_DS_temp**1.5))/(122.4+12.9*avg_molecular_weight+local_DS_temp))*0.0001
        K_DS=((9.379+0.016*avg_molecular_weight)*((1.8*local_DS_temp)**1.5))/(209.2+19.26*avg_molecular_weight+1.8*local_DS_temp)*0.0001
        K_AN=((9.379+0.016*avg_molecular_weight)*((1.8*local_AN_temp)**1.5))/(209.2+19.26*avg_molecular_weight+1.8*local_AN_temp)*0.0001

        #X_DS=2.57+(1914.5/local_DS_temp)+0.0095*avg_molecular_weight
        X_DS=3.448+(986.4/(local_DS_temp*1.8))+0.01*avg_molecular_weight
        X_AN=3.448+(986.4/(local_AN_temp*1.8))+0.01*avg_molecular_weight

        #Y_DS=1.11+0.004*X_DS
        Y_DS=2.447-0.2*X_DS
        Y_AN=2.447-0.2*X_AN
        
        exp_term_DS=X_DS*(local_DS_density**Y_DS)
        exp_term_AN=X_AN*(local_AN_density**Y_AN)


        local_gas_viscosity_DS=K_DS*math.exp(exp_term_DS)   
        local_gas_viscosity_AN=K_AN*math.exp(exp_term_AN)

        self.viscosities['DS viscosities'].append(local_gas_viscosity_DS)   
        self.viscosities['AN viscosities'].append(local_gas_viscosity_AN)   

    def calc_z_factors(self):
        for i in range(self.num_cells):
            self.calc_local_z_factor(i)
    
    def calc_volumes(self):
        for i in range(self.num_cells):
            self.calc_local_volume_per_mole(i)

    def calc_densities(self):
        for i in range(self.num_cells):
            self.calc_local_gas_density(i)
         
    def calc_viscosities(self):
        for i in range(self.num_cells):
            self.calc_local_gas_viscosity(i)

