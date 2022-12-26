from temperature import Temperature
from topology import Topology
from flow import Flow
import csv
import pandas as pd
import matplotlib as plt

def main():
    run_time=600


    Simulation=Flow()
    temp_over_time={0:Simulation.Ar}
    #velocity is const in our set up, only changes when flow rate does, which is not the case
    #static pressure is also const
    static_pressure={0:Simulation.static_pressure}
    liquid_sup_velocities={0:Simulation.liquid_sup_velocities}
    gas_sup_velocities={0:Simulation.gas_sup_velocities}
    z_factors={0:Simulation.z_factors}
    gas_densities={0:Simulation.gas_densities}
    gas_viscosities={0:Simulation.viscosities}
    reynolds_numbers={0:Simulation.reynolds_numbers}
    flow_types={0:Simulation.flow_types}
    slip_velocities={0:Simulation.slip_velocities}
    void_fractions={0:Simulation.void_fractions}
    friction_factors={0:Simulation.friction_factors}
    pressure_loss={0:Simulation.pressure_drops}

    for time in range(1,run_time):
        Simulation.calc_temp()
        temp_over_time[time]=Simulation.Ar

        Simulation.calc_sup_velocities()
        liquid_sup_velocities[time]=Simulation.liquid_sup_velocities
        gas_sup_velocities[time]=Simulation.gas_sup_velocities

        Simulation.calc_z_factors()
        z_factors[time]=Simulation.z_factors

        Simulation.calc_densities()
        gas_densities[time]=Simulation.gas_densities

        Simulation.calc_viscosities()
        gas_viscosities[time]=Simulation.viscosities

        Simulation.calc_reynolds_numbers()
        reynolds_numbers[time]=Simulation.reynolds_numbers

        Simulation.calc_flow_type()
        flow_types[time]=Simulation.flow_types
        slip_velocities[time]=Simulation.slip_velocities

        Simulation.calc_void_fraction()
        void_fractions[time]=Simulation.void_fractions

        Simulation.calc_friction_factor()
        friction_factors[time]=Simulation.friction_factors

        Simulation.calc_pressure_loss()
        pressure_loss[time]=Simulation.pressure_drops


if __name__ == "__main__":
    main()