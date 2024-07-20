#include <iostream>
#include <vector>
#include "GFUM.cpp"
#include "FUM.cpp"
#include "DFUM.cpp"
#include <cmath>
#include <fstream>


template <typename T>
void pvec(T &vec){
    for (size_t i = 0; i < vec.size(); i++)
    {
        std::cout << vec[i] << " ";
    }
    
}

void write_vector_to_file(const std::vector<double>& vector, const std::string& filename) {
    std::ofstream file;
    file.open(filename, std::ios_base::app); 
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    for (const double& value : vector) {
        file << value << ' ';
    }
    file << '\n';
    file.close();
}

std::vector<double> generateNonUniformArray(int size, double min_val, double max_val) {
    std::vector<double> result;
    result.reserve(size);

    double base = 2.0; // Adjust this base to change the spacing

    for (int i = 1; i <= size; ++i) {  // Start from 1 and go to size
        double exponent = static_cast<double>(i) / (size + 1);  // +1 to exclude 0 and 1
        double value = min_val + (max_val - min_val) * (std::pow(base, exponent) - 1) / (base - 1);
        result.push_back(value);
    }

    return result;
}


int main() {
    //double mass = 1.0;
    double v0ball = 0.5;
    double a0 = 1e-2;
    double phi0 = 0.0;

    std::cout << "1e-5, Dfum, diff" << '\n';
    //auto fermi = FUM2(1.0,1e-2, 1e-5,0);
    std::vector<double> mass = {1.0,1.0};
    std::vector<double> v0arr = {-0.1,-0.1};
    std::vector<double> x0arr = {0.2,1.0};
    auto fermi = GeneralFUM(mass,x0arr, v0arr, 1e-5,0);
    //auto fermi = DoubleFUM(1.0,1.0, -0.1,-0.1, 1e-2,0);
    std::vector<double> diffarr = fermi.diffusion(10,10,100);
    // // //write_vector_to_file(diffarr,"./data/FUM2_diff_01k_01k_1k_noabs.txt");
    pvec<std::vector<double>> (diffarr);
    int size = 32;
    double min_val = 0.0;
    double max_val = 1.0;

    

    
    

    if (0){
        std::vector<double> phi;
        std::vector<double> speed;

        auto v0arr = linspace(0.01001,0.4,22);
        auto phases = linspace(0,2*M_PI,22);

        for (size_t i = 0; i < v0arr.size(); ++i) {
            double v0 = v0arr[i];
            for (size_t k = 0; k < phases.size(); ++k) {
                double phase = phases[k];
                FUM2 fermi(1, -v0, a0, phase);
                //DoubleFUM fermi(1.0,1.0, -v0,-v0,  a0, phase);
        
                for (int j = 0; j < 160; ++j) {
                    bool hit = fermi.rwallhit();
                    if (hit==false) {
                        std::cout << v0 << ", " << phase << std::endl;
                        break;
                    }
                    phi.push_back(fermi.mw.phi);
                    speed.push_back(std::abs(fermi.ball.speed));
                }
            }
            std::cout << "done: " << static_cast<double>(i + 1) / v0arr.size() << std::endl;
        }

        // Write data to a text file
        std::ofstream outFile("./data/FUM2_ps_1e-2_22_22_160_2.txt");
        for (size_t i = 0; i < phi.size(); ++i) {
            outFile << std::fixed << std::setprecision(6) << phi[i] << " " << speed[i] << "\n";
        }
        outFile.close();
    }
    if (0){
        
        //FUM2 fermi(1, -0.00147272727272727, 1e-3, 1.92536051417124);
        //FUM fermi(1, 0.000165656565656566, 1e-2,0.698131700797732);
        std::vector<double> v0arr = {-0.1,-0.1};
        double a0=1e-2;
        //std::vector<double> x0arr = generateNonUniformArray(2,0.001, 0.89);
        std::vector<double> x0arr = {0.19,0.77};
        std::vector<double> mass = {0.01,1};
        auto fermi = DoubleFUM(0.01,1.0, -0.1,-0.1,a0,0);
        //auto fermi = GeneralFUM(mass, x0arr,v0arr,a0,0);
        std::setprecision(6);
        std::cout << fermi.mw.position << '\n';
        std::cout << fermi.ballr.position << '\n';
        std::cout << fermi.ballr.speed << '\n';
        std::cout << fermi.mw.speed << '\n';
        
        for (size_t i = 0; i < 10; i++)
        {
            std::cout << i << '\n';
            auto ch = fermi.collision();
            std::cout << ch << '\n';
            std::cout << ch << '\n';
            std::cout << fermi.mw.position << '\n';
            std::cout << fermi.ballr.position << '\n';
            std::cout << fermi.mw.speed << '\n';
            std::cout << fermi.ballr.speed << '\n';
            std::cout << fermi.balll.position << '\n';
            std::cout << "-----------------" << '\n';

        }
        
    }
    
    if (0){
        a0 = 1e-1;
        
        std::vector<double> v0arr = {-0.1,-0.1,-0.1,-0.1};
        std::vector<double> x0arr = generateNonUniformArray(8,0.001, 0.89);
        //std::vector<double> x0arr = {0.19,0.77};
        //std::cout << "rbhit, 1e-4, dfum, diff_1000" << '\n';
        int Niter = 100;
        int wdif = Niter/100;
        std::vector<double> mm = {1e-5,5e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.5};
        std::vector<double> darr = {};
        for (size_t i = 0; i < mm.size(); i++)
        {
            std::vector<double> mass = {mm[i],1.0,mm[i],1.0};
            //auto fermi = DoubleFUM(mm[i],1.0, -0.1,-0.1,a0,0);
            auto fermi = GeneralFUM(mass, x0arr,v0arr,a0,0);
            std::vector<double> diffarr = fermi.diffusion(100,100,Niter);
            darr.push_back(diffarr[wdif]);
        }
        std::vector<double> mass = {1,1,1,1};
        //auto fermi = DoubleFUM(1,1.0, -0.1,-0.1,a0,0);
        auto fermi = GeneralFUM(mass, x0arr,v0arr,a0,0);
        std::vector<double> diffarr = fermi.diffusion(100,100,Niter);
        darr.push_back(diffarr[wdif]);
        for (size_t i = mm.size(); i-- > 0;)
        {
            std::vector<double> mass = {1,mm[i],1,mm[i]};
            //auto fermi = DoubleFUM(1.0,mm[i], -0.1,-0.1,a0,0);
            auto fermi = GeneralFUM(mass, x0arr,v0arr,a0,0);
            std::vector<double> diffarr = fermi.diffusion(100,100,Niter);
            darr.push_back(diffarr[wdif]);
        }
        
        write_vector_to_file(darr,"./data/Gfum_01k_01k_Dm4_new2.txt");
    }
    return 0;
};