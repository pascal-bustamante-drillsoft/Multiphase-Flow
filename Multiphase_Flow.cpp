#define _USE_MATH_DEFINES
 
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

struct 
{
    int ir;
    int or;
    int oh;
    float well[2] = {3500,4500};
} Dimentions;


class Constants{
    public:
        float dt = 1;
        float dz = 100;
        float h1 = 10;
        float h2 = 10;
        float h3 = 10;
        float hef = 0.4;
        float l1 = 1;
        float l2 = 8.09;
        float l3 = 1;
        float c1 = 0.4;
        float q = 50;
        float GG = 0.0127;
        float supTemp = 60;
        float den = 10;
        float ir = 4;  //inner radus needs to be changable
        float Q1 = 0;   // UNKNOWN
        float height = 5000;
        
    
};

class Grid{
    public:
        float ir =  1;
        float or =  1;
        float oh =  1;
        float numCells = 50;
        vector<vector<float>>grid = {};  
        auto gridConst(float height, float dz);  
};

auto Grid::gridConst(float height, float dz){
    float z = 0;
    float ir = 4;
    float or = 5;
    vector<float> cell ={};
    while (numCells > 0){
        if (z == 3500){
            ir = 3.5;
            or = 5.5;
        }
        if (z == 4500){
            ir = 2.5;
            or = 6;
        }
        vector<float> cell = {z,ir, or};
        Grid::grid.push_back(cell);
    }
    return Grid::grid;
};


class Coefficients: private Constants{
    public:
        float A = (den*c1)/dt;
        float B = (den*c1*q)/(M_PI*(ir*ir)*dz);   //ir needs to change
        float C = (2*h1)/ir;
        float E = 1/(A+B+C);
        float D = 0;             //Q1 unknown
};

class Temperature: private Coefficients{
    vector<vector<float>> T0, T00;
    public:
        vector<vector<float>> temp;
        auto calcTemp(vector<vector<float>> Ar);
};

auto Temperature::calcTemp(vector<vector<float>> Ar){
    int length = sizeof(Ar[0]);
    float t1 = Ar[1][0];
    float t2 = Ar[0][0];
    float t3 = Ar[0][0];

    float t4 = ((t1*Coefficients::A)+(t2*Coefficients::B)+(t3*Coefficients::C)+Coefficients::D)*Coefficients::E;

    vector<float> res = {{t4}};

    for (int step = 1; step < length; step++){
        float t1 = Ar[1][step];
        float t2 = res[step-1];
        float t3 = Ar[0][step];
        
        float t4 = ((t1*Coefficients::A)+(t2*Coefficients::B)+(t3*Coefficients::C)+Coefficients::D)*Coefficients::E;
        res.push_back(t4);

    };
    Ar.erase(Ar.begin());
    Ar. push_back(res);
    return Ar;
};

int main(){
    vector<vector<float>> Ar = {{60,61.27,62.54,63.81,65.08,66.35,67.62},{60,61.27,62.54,63.81,65.08,66.35,67.62}};
    Temperature P;
    std::ofstream myfile;
    myfile.open("temperatures.csv");

    if (!myfile) {
        std::cerr << "can't open output file" << std::endl;
    }
    int max_time = 60;
    int j,i=0; 
    for(j=0; j < max_time; j++){
        Ar = P.calcTemp(Ar);
        myfile << Ar[0][0];
        for(i = 1; i < 7; i++){
            myfile <<", " << Ar[0][i];
        }
        myfile << "\n";
    }
    

    
    myfile.close();
    return 0;
}