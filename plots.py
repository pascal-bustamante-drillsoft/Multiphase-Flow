import glob 
import csv
import math
from iapws import IAPWS97
import matplotlib.pyplot as plt 
import numpy as np

class Plot():
    def __init__(self,paths={'path_temp_DS':'plots\temperatures2_DS.csv','path_temp_AN':'plots\temperatures2_AN.csv','path_temp_DS_wall':'plots\temperatures2_DS_wall.csv'}) -> None:
        """
        Get paths; Note: keep paths global (not relative)
        It could create a bug with PyInstaller
        """
        self.path_temp_DS=paths['path_temp_DS']
        self.path_temp_AN=paths['path_temp_AN']
        self.path_temp_DS_wall=paths['path_temp_DS_wall']
    
    def get_csv_data(self,file):
        """
        simple csv reader
        """
        table=[]
        with open(file,'r') as f:
            csv_reader=csv.reader(f)
            for row in csv_reader:
                for ele in row:
                    point=float(ele)
                    table.append(point)
        return table
    
    def plot_something(data, ax=None, **kwargs):
        ax = ax or plt.gca()
        # Do some cool data transformations...
        
        return ax.plot(data, **kwargs)



path_temp_DS='plots\temperatures2_DS.csv'
path_temp_AN='plots\temperatures2_AN.csv'
path_temp_DS_wall='plots\temperatures2_DS_wall.csv'

def get_temp(file):
    temperatures = []
    with open(file, "r") as f:
        csvreader=csv.reader(f)
        for row in csvreader:
            for x in row:
                temp = float(x)
                temperatures.append(temp)
    return temperatures








x = [n*100 for n in range(50)]

with open('temperatures2_DS.csv', mode='r') as f:
    csvreader=csv.reader(f)
    rows_ds=[]
    for row in csvreader:
        r=[]
        for z in row:
            temp = float(z)
            r.append(temp)
        rows_ds.append(r)

with open('temperatures2_AN.csv', mode='r') as f:
    csvreader=csv.reader(f)
    rows_an=[]
    for row in csvreader:
        r=[]
        for z in row:
            temp = float(z)
            r.append(temp)
        rows_an.append(r)

with open('temperatures2_AN.csv', mode='r') as f:
    csvreader=csv.reader(f)
    rows_an=[]
    for row in csvreader:
        r=[]
        for z in row:
            temp = float(z)
            r.append(temp)
        rows_an.append(r)
with open('temperatures2_DS_wall.csv', mode='r') as f:
    csvreader=csv.reader(f)
    rows_wall=[]
    for row in csvreader:
        r=[]
        for z in row:
            temp = float(z)
            r.append(temp)
        rows_wall.append(r)


y=rows_an[200]
z=rows_ds[200]
k=rows_wall[200]

plt.plot(y,x,label='annulus')  # Plot the chart
plt.plot(z,x,label='drill string')
plt.plot(k,x,label='drill string wall')
plt.legend()
plt.ylabel('temperature (F)')
plt.gca().invert_yaxis()
plt.xlabel('depth (ft)')
plt.title('Temperature')
plt.show()



