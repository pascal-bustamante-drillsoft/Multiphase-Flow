from temperature import Temperature
from topology.topology import Topology
import csv
import pandas as pd
import matplotlib as plt

def main():
    """
    oder:
    -initialize the grid (pd.DataFrame) via a topology (data)
    -
    """
    #initialize the topology first, then calc temp
    top=Topology()
    top.populate_radii()
    top.calc_cross_sections

    T = Temperature()
    Arr = T.Ar
    T.Ar = T.fahrenheit_to_kelvin(T.Ar)
    
    with open('temperatures2_DS.csv', 'w', encoding='UTF8', newline='') as f1:
        with open('temperatures2_DS_wall.csv', 'w', encoding='UTF8', newline='') as f2:
            with open('temperatures2_AN.csv', 'w', encoding='UTF8', newline='') as f3:

                writer1 = csv.writer(f1)
                writer2 = csv.writer(f2)
                writer3 = csv.writer(f3)

                for i in range(600):
                    writer1.writerow(T.kelvin_to_fahrenheit(Arr[0]))
                    writer2.writerow(T.kelvin_to_fahrenheit(Arr[1]))
                    writer3.writerow(T.kelvin_to_fahrenheit(Arr[2]))
                    T.calc_temp()
                   # print("in main ",T.kelvin_to_fahrenheit(T.Ar[2][0]))

if __name__ == "__main__":
    main()
