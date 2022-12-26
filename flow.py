from materials import Gas
from velocity import Velocity 
import math

'''
Flow of calculation for the flow:
    -
'''

class Flow(Gas,Velocity):
    def __init__(self, in_temp=60, dt=1, DS_flow_rate=0.031) -> None:
        super().__init__(in_temp=in_temp, dt=dt, DS_flow_rate=DS_flow_rate)

        self.flow_types={'DS flow types':[],'AN flow types':[]}
        self.slip_velocities={'DS slip velocities':[],'AN slip velocities':[]}
        self.void_fractions={'DS void fractions':[],'AN void fractions':[]}
        self.water_viscosity={'DS viscosity':[], 'AN viscosity':[]}
        self.reynolds_numbers={'DS reynolds numbers':[],'AN reynolds numbers':[]}
        self.pressure_drops={'DS pressure drops':[],'AN pressure drops':[]}
        self.friction_factors={'DS friction factors':[],'AN friction factors':[]}

        self.calc_reynolds_numbers()
        #print(self.reynolds_numbers)
        #print(self.liquid_sup_velocities)

        self.calc_flow_type()
        #print(self.flow_types)
        self.calc_void_fraction()
        #print(self.slip_velocities)
        self.calc_friction_factor()
        self.calc_pressure_loss()
        #print(self.pressure_drops)



    def calc_oil_viscosity(self,cell):  #not implemented there are constant values in Constants class
        """
        2 types of viscosuty:
        -one for saturated crude oil
        -the other for undersaturated crude oil
        """
        pass


    def calc_water_viscosity(self,cell):
        temp_DS=self.kelvin_to_fahrenheit(self.t_DS[cell])
        temp_AN=self.kelvin_to_fahrenheit(self.t_AN[cell])
        water_viscosity_DS=math.exp(1.003-(0.01479*temp_DS)+(0.00001982*(temp_DS**2)))
        water_viscosity_AN=math.exp(1.003-(0.01479*temp_AN)+(0.00001982*(temp_AN**2)))

        if len(self.water_viscosity['DS viscosity']) != cell:
            raise Exception('Indexing problem in DS viscosity calculation')
        if len(self.water_viscosity['AN viscosity']) != cell:
            raise Exception('Indexing problem in AN viscosity calculation')

        self.water_viscosity['DS viscosity'].append(water_viscosity_DS)
        self.water_viscosity['AN viscosity'].append(water_viscosity_AN)

    def calc_local_reynolds_numbers(self,cell):
        '''
        It normaly has several variables like:
        [density, viscosity,characteristic length]
        but for us the only variable is characteristic length
        '''
        characteristic_length_DS=2*0.0254*(self.inner_radii[cell])
        characteristic_length_AN=2*0.0254*(self.bore_radii[cell]-self.outer_radii[cell])

        velocity_DS=self.liquid_sup_velocities['DS superficial velocities'][cell]
        velocity_AN=self.liquid_sup_velocities['AN superficial velocities'][cell]

        reynolds_number_DS=self.den_DS*velocity_DS*characteristic_length_DS/self.mat_viscosity
        reynolds_number_AN=self.den_AN*velocity_AN*characteristic_length_AN/self.mat_viscosity

        #print(self.den_DS,velocity_DS,characteristic_length_DS,self.mat_viscosity)

        self.reynolds_numbers['DS reynolds numbers'].append(reynolds_number_DS)
        self.reynolds_numbers['AN reynolds numbers'].append(reynolds_number_AN)


    def calc_flow_type_DS(self,cell):
        '''
        This comes from Baojiang Sun
        Multiphase Flow in Oil and Gas Well Drilling, 2016
        pages 103-106
        '''
        v_sg_DS=self.gas_sup_velocities['DS superficial velocities'][cell]

        v_sl_DS=self.liquid_sup_velocities['DS superficial velocities'][cell]

        gas_den_DS=self.gas_densities['DS densities'][cell]
        
        d_t=self.inches_to_meters(self.inner_radii[cell])
        d_c=self.inches_to_meters(self.outer_radii[cell])
        r=d_t/d_c  

        '''
        Correction coefs not known
        '''
        k1=1
        k2=1

        '''
        The next become the bubble drift velocity 
        depending on the flow type
        '''
       # print(1.5(self.grav*self.surface_tension*(self.den_DS-gas_den_DS)/(self.den_DS**2))**0.25)
        v_bubble_bubbly_DS=1.5*((self.grav*self.surface_tension*(self.den_DS-gas_den_DS)/(self.den_DS**2))**0.25)
        v_taylor_bubble_slug_DS=(0.3+0.22*(r))*((self.grav*2*(d_c-d_t)*(self.den_DS-gas_den_DS)/self.den_DS)**0.5)
        v_taylor_bubble_churn_DS=(0.3+0.22*(r))*((self.grav*2*(d_c-d_t)*(self.den_AN-gas_den_DS)/self.den_DS)**0.5)

        '''
        Define conditions for types of flow
        then we asign the corresponding slig velicities
        '''
        bubbly_cond_DS=k1*(0.429*v_sl_DS+0.357*v_bubble_bubbly_DS)     #This might be wrong, there might be a typo in textbook page 103

        slug_cond_DS=0.429*v_sl_DS+0.357*v_taylor_bubble_slug_DS

        churn_cond_1_1_DS=25.4*math.log(self.den_DS*(v_sl_DS**2))-38.9
        churn_cond_1_2_DS=74.4
        churn_cond_2_1_DS=0.0051*((self.den_DS*(v_sl_DS**2))**1.7)
        churn_cond_2_2_DS=14.4
        churn_cond_3_DS=k2*((self.surface_tension*self.grav*((self.den_DS-gas_den_DS)**0.333)/(gas_den_DS**2))**0.25)

        annular_cond_DS=k2*((self.surface_tension*self.grav*((self.den_DS-gas_den_DS)**0.333)/(gas_den_DS**2))**0.25)
        #print(slug_cond_DS, v_taylor_bubble_slug_DS)

        if v_sg_DS<bubbly_cond_DS:
            flow_DS='bubbly'
            local_slip_velocity_DS=v_bubble_bubbly_DS
            self.flow_types['DS flow types'].append(flow_DS)
            self.slip_velocities['DS slip velocities'].append(local_slip_velocity_DS)

        elif v_sg_DS>slug_cond_DS:
            flow_DS='slug'
            local_slip_velocity_DS=v_taylor_bubble_slug_DS
            self.flow_types['DS flow types'].append(flow_DS)
            self.slip_velocities['AN slip velocities'].append(local_slip_velocity_DS)

        elif gas_den_DS*(v_sg_DS**2)>churn_cond_1_1_DS and self.den_DS*(v_sl_DS**2)>churn_cond_1_2_DS:

            if gas_den_DS*(v_sg_DS**2)>churn_cond_2_1_DS and self.den_DS*(v_sl_DS**2)>churn_cond_2_2_DS:

                if v_sg_DS<churn_cond_3_DS:
                    flow_DS='churn'
                    local_slip_velocity_DS=v_taylor_bubble_churn_DS
                    self.flow_types['DS flow types'].append(flow_DS)
                    self.slip_velocities['AN slip velocities'].append(local_slip_velocity_DS)
        
        elif v_sg_DS>annular_cond_DS:
            flow_DS='annular'
            local_slip_velocity_DS=v_sg_DS
            self.flow_types['DS flow types'].append(flow_DS)
            self.slip_velocities['AN slip velocities'].append(local_slip_velocity_DS)

        else:
            raise Exception('flow type error, none of the flow type condition were met')

    def calc_flow_type_AN(self,cell):
        '''
        This comes from Baojiang Sun
        Multiphase Flow in Oil and Gas Well Drilling, 2016
        pages 103-106
        '''
        v_sg_AN=self.gas_sup_velocities['AN superficial velocities'][cell]

        v_sl_AN=self.liquid_sup_velocities['AN superficial velocities'][cell]

        gas_den_AN=self.gas_densities['AN densities'][cell]
        
        d_c=self.bore_radii[cell]
        d_t=self.outer_radii[cell]
        r=d_t/d_c  

        '''
        Correction coefs not known
        '''
        k1=5
        k2=1

        '''
        The next become the bubble drift velocity 
        depending on the flow type
        '''
        v_bubble_bubbly_AN=1.5*((self.grav*self.surface_tension*(self.den_AN-gas_den_AN)/(self.den_AN**2))**0.25)

        v_taylor_bubble_slug_AN=(0.3+0.22*(r))*((self.grav*2*(d_c-d_t)*(self.den_AN-gas_den_AN)/self.den_AN)**0.5)

        v_taylor_bubble_churn_AN=(0.3+0.22*(r))*((self.grav*2*(d_c-d_t)*(self.den_AN-gas_den_AN)/self.den_AN)**0.5)

        '''
        Define conditions for types of flow
        then we asign the corresponding slig velicities
        '''
        bubbly_cond_AN=k1*(0.429*v_sl_AN+0.357*v_bubble_bubbly_AN)     #This might be wrong, corrected typo in textbook page 103

        slug_cond_AN=0.429*v_sl_AN+0.357*v_taylor_bubble_slug_AN

        churn_cond_1_1_AN=25.4*math.log(self.den_AN*(v_sl_AN**2))-38.9
        churn_cond_1_2_AN=74.4
        churn_cond_2_1_AN=0.0051*((self.den_AN*(v_sl_AN**2))**1.7)
        churn_cond_2_2_AN=14.4
        churn_cond_3_AN=k2*((self.surface_tension*self.grav*((self.den_AN-gas_den_AN)**0.333)/(gas_den_AN**2))**0.25)

        annular_cond_AN=k2*((self.surface_tension*self.grav*((self.den_AN-gas_den_AN)**0.333)/(gas_den_AN**2))**0.25)
        #print(v_sg_AN,bubbly_cond_AN,slug_cond_AN,annular_cond_AN)

        if v_sg_AN<bubbly_cond_AN:
            flow_AN='bubbly'
            local_slip_velocity_AN=v_bubble_bubbly_AN
            self.flow_types['AN flow types'].append(flow_AN)
            self.slip_velocities['AN slip velocities'].append(local_slip_velocity_AN)

        elif v_sg_AN>slug_cond_AN:
            flow_AN='slug'
            local_slip_velocity_AN=v_taylor_bubble_slug_AN
            self.flow_types['AN flow types'].append(flow_AN)
            self.slip_velocities['AN slip velocities'].append(local_slip_velocity_AN)

        elif gas_den_AN*(v_sg_AN**2)>churn_cond_1_1_AN and self.den_AN*(v_sl_AN**2)>churn_cond_1_2_AN:

            if gas_den_AN*(v_sg_AN**2)>churn_cond_2_1_AN and self.den_AN*(v_sl_AN**2)>churn_cond_2_2_AN:

                if v_sg_AN<churn_cond_3_AN:
                    flow_AN='churn'
                    local_slip_velocity_AN=v_taylor_bubble_churn_AN
                    self.flow_types['AN flow types'].append(flow_AN)
                    self.slip_velocities['AN slip velocities'].append(local_slip_velocity_AN)
        
        elif v_sg_AN>annular_cond_AN:
            flow_AN='annular'
            local_slip_velocity_AN=v_sg_AN
            self.flow_types['AN flow types'].append(flow_AN)
            self.slip_velocities['AN slip velocities'].append(local_slip_velocity_AN)

        else:
            raise Exception('flow type calculation error, none of the flow type condition were met')

    
    def calc_local_void_fraction_DS(self,cell):
        '''
        Using P1V1Z1T0=P0V0Z0T1
        surface tension will be constant as 30 dyn/cm
        '''
        surface_tension=self.surface_tension
        d_t=self.inches_to_meters(self.outer_radii[cell])
        d_c=self.inches_to_meters(self.bore_radii[cell])
        r=d_t/d_c                                                           #was taken out of DS calculations

        gas_den_DS=self.gas_densities['DS densities'][cell]
        gas_viscosity_DS=self.viscosities['DS viscosities'][cell]

        flow=self.flow_types['DS flow types'][cell]

        v_sg_DS=self.gas_sup_velocities['DS superficial velocities'][cell]
        v_sl_DS=self.liquid_sup_velocities['DS superficial velocities'][cell]

        v_m_DS=v_sg_DS+v_sl_DS


        if flow == 'bubbly': 
            c_DS=1.2

            v_gas_bubble_DS=1.5*((self.grav*surface_tension*(self.den_DS-gas_den_DS)/(self.den_DS**2))**0.25)
            local_void_fraction_DS=(v_sg_DS)/(c_DS*v_m_DS+v_gas_bubble_DS)
            #print(local_void_fraction_DS)

            self.void_fractions['DS void fractions'].append(local_void_fraction_DS)
        
        elif flow == 'slug': 
            
            c_DS=1.18       #got rid of r term for cylinder (its a annulus term)
            
            v_gas_taylor_bubble_DS=(0.3)*((self.grav*2*(d_c-d_t)*(self.den_DS-gas_den_DS)/self.den_DS)**0.5)
            local_void_fraction_DS=(v_sg_DS)/(c_DS*v_m_DS+v_gas_taylor_bubble_DS)
            #print(local_void_fraction_DS)


            self.void_fractions['DS void fractions'].append(local_void_fraction_DS)
        
        elif flow == 'churn':
            c_DS=1.15

            v_gas_taylor_bubble_DS=(0.3)*((self.grav*2*(d_c-d_t)*(self.den_AN-gas_den_DS)/self.den_DS)**0.5)
            local_void_fraction_DS=(v_sg_DS)/(c_DS*v_m_DS+v_gas_taylor_bubble_DS)
            #print(local_void_fraction_DS)


            self.void_fractions['DS void fractions'].append(local_void_fraction_DS)

        elif flow == 'annular':
            v_critical_DS=(v_sg_DS*gas_viscosity_DS*((gas_den_DS/self.den_DS)**0.5))/surface_tension
            if (v_critical_DS*10000)<=4:
                E=0.0055*((v_critical_DS*10000)**2.86)
            else:
                E=0.857*math.log(v_critical_DS*10000)-0.2

            local_void_fraction_DS=v_sg_DS/(v_sg_DS+E*v_sl_DS)
            #print(local_void_fraction_DS)

            
            self.void_fractions['DS void fractions'].append(local_void_fraction_DS)

        else:
            raise Exception ("Void fraction error, no flow valid flow type")
    
    def calc_local_void_fraction_AN(self,cell):
        '''
        surface tension will be constant as 30 dyn/cm
        r=ratio of tubular diameter to  
        '''
        surface_tension=self.surface_tension
        d_t=self.inches_to_meters(self.outer_radii[cell])
        d_c=self.inches_to_meters(self.bore_radii[cell])
        r=d_t/d_c  

        gas_den_AN=self.gas_densities['AN densities'][cell]

        gas_viscosity_AN=self.viscosities['AN viscosities'][cell]

        v_sg_AN=self.gas_sup_velocities['AN superficial velocities'][cell]
        v_sl_AN=self.liquid_sup_velocities['AN superficial velocities'][cell]

        v_m_AN=v_sg_AN+v_sl_AN

        flow=self.flow_types['DS flow types'][cell]

        if flow == 'bubbly': 
            c_AN=1.2+0.371*r

            v_gas_bubble_AN=1.5*((self.grav*surface_tension*(self.den_AN-gas_den_AN)/(self.den_AN**2))**0.25)
            local_void_fraction_AN=(v_sg_AN)/(c_AN*v_m_AN+v_gas_bubble_AN)

            self.void_fractions['AN void fractions'].append(local_void_fraction_AN)
        
        elif flow == 'slug': 
            c_AN=1.18+0.9*r
            
            v_gas_taylor_bubble_AN=(0.3+0.22*(r))*((self.grav*2*(d_c-d_t)*(self.den_AN-gas_den_AN)/self.den_AN)**0.5)
            local_void_fraction_AN=(v_sg_AN)/(c_AN*v_m_AN+v_gas_taylor_bubble_AN)

            self.void_fractions['AN void fractions'].append(local_void_fraction_AN)
        
        elif flow == 'churn':
            c_AN=1.15

            v_gas_taylor_bubble_AN=(0.3+0.22*(r))*((self.grav*2*(d_c-d_t)*(self.den_AN-gas_den_AN)/self.den_AN)**0.5)
            local_void_fraction_AN=(v_sg_AN)/(c_AN*v_m_AN+v_gas_taylor_bubble_AN)

            self.void_fractions['AN void fractions'].append(local_void_fraction_AN)

        elif flow == 'annular':
            v_critical_AN=(v_sg_AN*gas_viscosity_AN*((gas_den_AN/self.den_AN)**0.5))/surface_tension

            if (v_critical_AN*10000)<=4:
                E=0.0055*((v_critical_AN*10000)**2.86)
            else:
                E=0.857*math.log(v_critical_AN*10000)-0.2
                
            local_void_fraction_AN=v_sg_AN/(v_sg_AN+E*v_sl_AN) 
            self.void_fractions['AN void fractions'].append(local_void_fraction_AN) 

        else:
            raise Exception ("Void fraction error, no flow valid flow type")

        '''
        The rest is based on
        Hassan and Kabir's paper
        Performance of a two-phase gas/liquid flow model in vertical wells
        1990
        '''

    def calc_local_friction_factor_DS(self,cell):
        '''
        Based on An Explicit Equation for Friction Factors in Pipes. Industrial & Engineering Chemistry Fundamentals
        Chen (1979)
        This is for a cylinder
        '''
        D_DS=2*self.inner_radii[cell]
        Re_DS=self.reynolds_numbers['DS reynolds numbers'][cell]

        local_relative_pipe_roughness=self.pipe_roughness/D_DS
        A=((local_relative_pipe_roughness**1.1098)/2.8257)+((7.149/Re_DS)** 0.8981)
        C=5.0452*math.log(A)/Re_DS
        D=(local_relative_pipe_roughness/3.7065)
        E=abs(D-C)
        B=math.log(E)
        f=1/((4*B)**2)

        self.friction_factors['DS friction factors'].append(f)

    def calc_local_friction_factor_AN(self,cell):
        '''
        Based on An Explicit Equation for Friction Factors in Pipes. Industrial & Engineering Chemistry Fundamentals
        Chen (1979)
        This is for an annulus
        '''
        D_AN=2*(self.bore_radii[cell]+self.outer_radii[cell])       #adding peripheries of annulus
        Re_AN=self.reynolds_numbers['AN reynolds numbers'][cell]

        local_relative_pipe_roughness=self.pipe_roughness/D_AN
        A=((local_relative_pipe_roughness**1.1098)/2.8257)+((7.149/Re_AN)** 0.8981)
        C=5.0452*math.log(A)/Re_AN
        D=(local_relative_pipe_roughness/3.7065)
        E=abs(D-C)
        B=math.log(E)
        f=1/((4*B)**2)

        self.friction_factors['AN friction factors'].append(f)

    def calc_local_pressure_loss_DS(self,cell):
        '''
        Based on Hassan and Kabir's paper
        Performance of a two-phase gas/liquid flow model in vertical wells
        1990
        '''
        flow=self.flow_types['DS flow types'][cell]
        f=self.friction_factors['DS friction factors'][cell]
        v_sg_DS=self.gas_sup_velocities['DS superficial velocities'][cell]
        v_sl_DS=self.liquid_sup_velocities['DS superficial velocities'][cell]
        gas_den_DS=self.gas_densities['DS densities'][cell]
        void_fraction=self.void_fractions['DS void fractions'][cell]
        Re_DS=self.reynolds_numbers['DS reynolds numbers'][cell]
        gas_viscosity_DS=self.viscosities['DS viscosities'][cell]
        D=2*self.inner_radii[cell]

        v_m_DS=v_sg_DS+v_sl_DS
        mat_den_DS=gas_den_DS*void_fraction+self.den_DS*(1-void_fraction)


        if flow == 'bubbly':
            local_frictional_pressure_drop=1000*2*f*(v_m_DS**2)*mat_den_DS/(D)                #g_c in paper =1 in SI units
            if cell == 0:
                self.pressure_drops['DS pressure drops'].append(local_frictional_pressure_drop)
            else:
                #print(self.pressure_drops)
                local_frictional_pressure_drop+=self.pressure_drops['DS pressure drops'][cell-1]
                self.pressure_drops['DS pressure drops'].append(local_frictional_pressure_drop)

        elif flow == 'slug' or flow == 'churn':
            local_frictional_pressure_drop=1000*2*f*(v_m_DS**2)*self.den_DS*(1-void_fraction)/(D)               #g_c in paper =1 in SI units
            if cell == 0:
                self.pressure_drops['DS pressure drops'].append(local_frictional_pressure_drop)
            else:
                local_frictional_pressure_drop+=self.pressure_drops['DS pressure drops'][cell-1]
                self.pressure_drops['DS pressure drops'].append(local_frictional_pressure_drop)

        elif flow == 'annular':
            f=0.079*(1+75*(1-void_fraction))/(Re_DS**0.25)                              #friction factor for annular flow
            v_critical_DS=(v_sg_DS*gas_viscosity_DS*((gas_den_DS/self.den_DS)**0.5))/self.surface_tension

            if (v_critical_DS*10000)<=4:
                E=0.0055*((v_critical_DS*10000)**2.86)
            else:
                E=0.857*math.log(v_critical_DS*10000)-0.2
                
            core_den_DS=((v_sg_DS*gas_den_DS)+(E*self.den_DS*v_sl_DS))/(v_sg_DS+E*v_sg_DS)      #there might be something wrong with the denominator (double v_sg)
            local_frictional_pressure_drop=1000*2*f*core_den_DS*((v_sg_DS/void_fraction)**2)/(D)     #g_c in paper =1 in SI units

            if cell == 0:
                self.pressure_drops['DS pressure drops'].append(local_frictional_pressure_drop)
            else:
                local_frictional_pressure_drop+=self.pressure_drops['DS pressure drops'][cell-1]
                self.pressure_drops['DS pressure drops'].append(local_frictional_pressure_drop)

        else:
            raise Exception('Frictional pressure loss calculation error, flow type must be calculated prior') 

    def calc_local_pressure_loss_AN(self,cell):
        '''
        Based on Hassan and Kabir's paper
        Performance of a two-phase gas/liquid flow model in vertical wells
        1990
        '''
        flow=self.flow_types['AN flow types'][cell]
        f=self.friction_factors['AN friction factors'][cell]
        v_sg_AN=self.gas_sup_velocities['AN superficial velocities'][cell]
        v_sl_AN=self.liquid_sup_velocities['AN superficial velocities'][cell]
        gas_den_AN=self.gas_densities['AN densities'][cell]
        void_fraction=self.void_fractions['AN void fractions'][cell]
        Re_AN=self.reynolds_numbers['AN reynolds numbers'][cell]
        gas_viscosity_AN=self.viscosities['AN viscosities'][cell]
        D=2*self.inner_radii[cell]

        v_m_AN=v_sg_AN+v_sl_AN
        mat_den_AN=gas_den_AN*void_fraction+self.den_AN*(1-void_fraction)


        if flow == 'bubbly':
            local_frictional_pressure_drop=1000*2*f*(v_m_AN**2)*mat_den_AN/(D)                #g_c in paper =1 in SI units
            if cell == 0:
                self.pressure_drops['AN pressure drops'].append(local_frictional_pressure_drop)
            else:
                local_frictional_pressure_drop+=self.pressure_drops['AN pressure drops'][cell-1]
                self.pressure_drops['AN pressure drops'].append(local_frictional_pressure_drop)

        elif flow == 'slug' or flow == 'churn':
            local_frictional_pressure_drop=1000*2*f*(v_m_AN**2)*self.den_AN(1-void_fraction)/(D)               #g_c in paper =1 in SI units
            if cell == 0:
                self.pressure_drops['AN pressure drops'].append(local_frictional_pressure_drop)
            else:
                local_frictional_pressure_drop+=self.pressure_drops['AN pressure drops'][cell-1]
                self.pressure_drops['AN pressure drops'].append(local_frictional_pressure_drop)

        elif flow == 'annular':
            f=0.079*(1+75*(1-void_fraction))/(Re_AN**0.25)                              #friction factor for annular flow
            v_critical_AN=(v_sg_AN*gas_viscosity_AN*((gas_den_AN/self.den_AN)**0.5))/self.surface_tension

            if (v_critical_AN*10000)<=4:
                E=0.0055*((v_critical_AN*10000)**2.86)
            else:
                E=0.857*math.log(v_critical_AN*10000)-0.2
                
            core_den_AN=((v_sg_AN*gas_den_AN)+(E*self.den_AN*v_sl_AN))/(v_sg_AN+E*v_sg_AN)      #there might be something wrong with the denominator (double v_sg)
            local_frictional_pressure_drop=1000*2*f*core_den_AN*((v_sg_AN/void_fraction)**2)/(D)     #g_c in paper =1 in SI units

            if cell == 0:
                self.pressure_drops['AN pressure drops'].append(local_frictional_pressure_drop)
            else:
                local_frictional_pressure_drop+=self.pressure_drops['AN pressure drops'][cell-1]
                self.pressure_drops['AN pressure drops'].append(local_frictional_pressure_drop)

        else:
            raise Exception('Frictional pressure loss calculation error, flow type must be calculated prior') 

    def calc_reynolds_numbers(self):
        for cell in range(self.num_cells):
            self.calc_local_reynolds_numbers(cell)

    def calc_flow_type(self):
        for cell in range(self.num_cells):
            self.calc_flow_type_DS(cell)
            self.calc_flow_type_AN(cell)

    def calc_void_fraction(self):
        for cell in range(self.num_cells):
            self.calc_local_void_fraction_DS(cell)
            self.calc_local_void_fraction_AN(cell)

    def calc_friction_factor(self):
        for cell in range(self.num_cells):
            self.calc_local_friction_factor_DS(cell)
            self.calc_local_friction_factor_AN(cell)

    def calc_pressure_loss(self):
        for cell in range(self.num_cells):
            self.calc_local_pressure_loss_DS(cell)
            self.calc_local_pressure_loss_AN(cell)
