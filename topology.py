import math

class Topology(): #add topology to constants instead
    """
    Creates the radial coordinates for the bore well. 
    It assumes they are cylinders with different radii. 
    Optional inputs are of the form of list of tuples [(depth,new radii)}]
    """
    def __init__(self, inner_radii_boundaries=[(0,4),(3500,3.5),(4500,2.5)], outer_radii_boundaries=[(0,5),(3500,5.5),(4500,6)], bore_wall_boundaries=[(0,10),(3500,8.5)], height=5000.0, dz=100): #each 
        inner_sorted=sorted(inner_radii_boundaries)
        outer_sorted=sorted(outer_radii_boundaries)
        bore_sorted=sorted(bore_wall_boundaries)

        #global properties
        self.borders=[inner_sorted,outer_sorted,bore_sorted]
        self.num_cylinders=len(self.borders)
        self.dz=dz
        self.height=height
        self.num_cells=int(height/dz)

        #edge case
        if dz == 0 or height == 0:
            raise Exception('Cell height or max depth cannot be 0') 
        if self.num_cylinders != 3:
            raise Exception('topological error in class Grid constructor')

        #local properties
        self.inner_radii=[]
        self.outer_radii=[]
        self.bore_radii=[]
        #self.radii=[self.inner_radii,self.outer_radii,self.bore_radii]
        self.DS_areas=[]
        self.AN_areas=[]
        self.pipe_areas=[]

        #dictionary to input into pd df
        self.topology={'inner radii': self.inner_radii,
                'outer radii': self.outer_radii,
                'bore radii': self.bore_radii,
                'DS areas': self.DS_areas,
                'AN areas': self.AN_areas,
                'pipe areas': self.pipe_areas
            }

    """ 
    Methods that populate the radii, 
    seperated by cylinder
    """
    def populate_inner_radii(self):
        inner_radii_borders=self.borders[0]

        for i in range(1,len(inner_radii_borders)):
            boundary_top=inner_radii_borders[i-1]
            boundary_bottom=inner_radii_borders[i]
            depth_start=boundary_top[0]
            depth_end=boundary_bottom[0]
            radius=boundary_top[1]

            depth_range=depth_end-depth_start
            num_cells=int(depth_range//100)
            temp_radii=[radius]*num_cells
            self.inner_radii=self.inner_radii+temp_radii
            
            if i==len(inner_radii_borders)-1:
                depth_start=boundary_bottom[0]
                depth_end=self.height
                radius=boundary_bottom[1]

                depth_range=depth_end-depth_start
                num_cells=int(depth_range//100)
                temp_radii=[radius]*num_cells
                self.inner_radii=self.inner_radii+temp_radii

        self.topology['inner radii']=self.inner_radii

    def populate_outer_radii(self):
        outer_radii_borders=self.borders[1]

        for i in range(1,len(outer_radii_borders)):
            boundary_top=outer_radii_borders[i-1]
            boundary_bottom=outer_radii_borders[i]
            depth_start=boundary_top[0]
            depth_end=boundary_bottom[0]
            radius=boundary_top[1]

            depth_range=depth_end-depth_start
            num_cells=int(depth_range//100)
            temp_radii=[radius]*num_cells
            self.outer_radii=self.outer_radii+temp_radii

            if i==len(outer_radii_borders)-1:
                depth_start=boundary_bottom[0]
                depth_end=self.height

                depth_range=depth_end-depth_start
                num_cells=int(depth_range//100)
                temp_radii=[radius]*num_cells
                self.outer_radii=self.outer_radii+temp_radii

        self.topology['outer radii']=self.outer_radii
    
    def populate_bore_radii(self):
        bore_radii_borders=self.borders[2]

        for i in range(1,len(bore_radii_borders)):
            boundary_top=bore_radii_borders[i-1]
            boundary_bottom=bore_radii_borders[i]
            depth_start=boundary_top[0]
            depth_end=boundary_bottom[0]
            radius=boundary_top[1]

            depth_range=depth_end-depth_start
            num_cells=int(depth_range//100)
            temp_radii=[radius]*num_cells
            self.bore_radii=self.bore_radii+temp_radii

            if i==len(bore_radii_borders)-1:
                depth_start=boundary_bottom[0]
                depth_end=self.height

                depth_range=depth_end-depth_start
                num_cells=int(depth_range//100)
                temp_radii=[radius]*num_cells
                self.bore_radii=self.bore_radii+temp_radii
        
        self.topology['bore radii']=self.bore_radii    

    """
    Methods that calculate and populate the cross sectional areas,
    seperated by cylinder
    """  
    def calc_DS_cross_sectional_area(self):
        for radius in self.inner_radii:
            area = math.pi*(radius**2)
            self.DS_areas.append(area)

        self.topology['DS areas']=self.DS_areas

    def calc_AN_cross_sectional_area(self):
        #print(len(self.bore_radii))
        for i in range(self.num_cells):   #indexing error here 
            outer_radius=self.outer_radii[i]
            bore_radius=self.bore_radii[i]

            area=math.pi*((bore_radius**2)-(outer_radius**2))
            self.AN_areas.append(area)
        
        self.topology['AN areas']=self.AN_areas

    def calc_pipe_cross_section(self):
        for i in range(self.num_cells):
            inner_radius=self.inner_radii[i]
            outer_radius=self.outer_radii[i]

            area=math.pi*((inner_radius**2)-(outer_radius**2))
            self.pipe_areas.append(area)

        self.topology['pipe areas']=self.pipe_areas

    def populate_radii(self):
        self.populate_inner_radii()
        self.populate_outer_radii()
        self.populate_bore_radii()

    def calc_cross_sections(self):
        self.calc_DS_cross_sectional_area()
        self.calc_AN_cross_sectional_area()
        self.calc_pipe_cross_section()

