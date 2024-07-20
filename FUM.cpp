
#include "helper.h"



class FUM {
public:
    FUM(double mass, double v0ball, double a0, double phi0) 
        : tnow(0.0), mw(a0, tnow, phi0), leftwall(0.0),ball(mass, v0ball,0) {
        double poss = mw.position;
        if (a0 > std::abs(v0ball) && (phi0 < M_PI/2 || phi0 > 3*M_PI/2)) {
            poss = mw.position-1e-6;
        }
        ball.update(v0ball, poss);
    }
    
    double tnow;
    MovingWall mw;
    Ball ball;
    double leftwall;

    bool rwcollision() {
        auto fmu = [this](double dtnew) {
            return mw.amp * std::cos(dtnew + mw.phi0 + tnow) - (ball.position - 1) + ball.speed * dtnew;
        };

        std::vector<double> roots = find_roots(1e-6, 2 * M_PI, 100, fmu);
        if (roots.empty()) {
            return false;
        } else {
            tnow = *std::min_element(roots.begin(), roots.end()) + tnow;
            mw.update(tnow);
            ball.update(std::abs(2 * mw.speed - ball.speed), mw.position);
            return true;
        }
    }

    bool lwcollision() {
        double T = (ball.position / ball.speed) + (1 - mw.amp) / ball.speed;
        auto fxi = [this, T](double dtnew) {
            return mw.amp * std::cos(dtnew + mw.phi0 + tnow + T) - ball.speed * dtnew + mw.amp;
        };

        std::vector<double> roots = find_roots(1e-15, 2 * M_PI, 100, fxi);
        if (roots.empty()) {
            return false;
        } else {
            tnow = *std::min_element(roots.begin(), roots.end()) + tnow + T;
            mw.update(tnow);
            ball.update(std::abs(2 * mw.speed - ball.speed), mw.position);
            return true;
        }
    }


    bool rwallhit() {
        bool rhw = rwcollision();
        if (rhw) {
            return true;
        } else {
            bool lcol = lwcollision();
            if (!lcol) {
                std::cerr << "No solution after, check why" << std::endl;
                return false;
            }
            return true;
        }
    }
    std::vector<double> diffusion(int Nv0, int Nph, int Niter)  {
        auto v0 = linspace(0.01*mw.amp, 0.02*mw.amp, Nv0); 
        auto phases = linspace(0.0, 2 * M_PI, Nph);
        int NP = Nv0 * Nph;
        std::vector<double> narr(100, 1.0);
        double diffu = 0;
        std::vector<double> intermediate_diffusions(100, 0.0);
        double dif = 0;

        for (size_t i = 0; i < v0.size(); i++) {
            for (size_t j = 0; j < phases.size(); j++) {
                FUM fermi(1.0, v0[i],mw.amp, phases[j]); 
                int k = 0;

                for (int iter = 1; iter <= Niter; iter++) {
                    int hit = fermi.rwallhit();
                    if (hit == 0) {
                        std::cout << v0[i] << ", " << phases[j] << std::endl;
                        NP -= 1;
                        break;
                    }

                    if (iter % 10 == 0 || iter==1) {
                        double speed_diff = std::pow(fermi.ball.speed- (v0[i]), 2);
                        narr[k] = static_cast<double>(iter);
                        intermediate_diffusions[k] += speed_diff;
                        k += 1;
                    }
                }
                diffu += std::pow(fermi.ball.speed - (v0[i]), 2);
            }
            std::cout << "done: " << static_cast<double>(i + 1) / v0.size() << std::endl;
        }
        narr[0]=1;
        for (size_t i = 0; i < intermediate_diffusions.size(); i++)
        {
           intermediate_diffusions[i] /= (narr[i]*2*NP);
        }
        std::cout << (diffu/(Niter*2*NP)) << '\n';
        return intermediate_diffusions;
    }

};

//SECOND IMPLEMENTATION

class FUM2 : public Observables {
public:
    FUM2(double mass, double v0ball, double a0, double phi0)
        : tnow(0.0), mw(a0, tnow, phi0), leftwall(0.0),ball(mass, v0ball,0) {
        double poss = mw.position;
        if (a0 > std::abs(v0ball)){ //&& (phi0 < M_PI / 2 || phi0 > 3 * M_PI / 2)) {
            poss = 1-1.1*a0;
        }
        poss = 1-1.1*a0;
        ball.update(v0ball, poss);
    }
    Ball ball;
    MovingWall mw;
    double tnow;
    double leftwall;
    
    double trwall() {
        // Case 1: ball.position >= (1 - mw.amp)
        if (ball.position >= (1 - mw.amp)) {
            auto fmu = [&](double dtnew) {
                return mw.amp * std::cos(dtnew + mw.phi0 + tnow)
                    - (ball.position - 1) - ball.speed * dtnew;
            };

            // Find roots using find_roots function
            std::vector<double> roots = find_roots(1e-6, 2 * M_PI, 100, fmu);

            // Determine tmin based on roots
            double tmin = std::numeric_limits<double>::infinity();
            for (double root : roots) {
                if (root > 0 && root < tmin) {
                    tmin = root;
                    
                }
            }
            
            return tmin;
        }
        // Case 2: ball.speed > 0
        else if (ball.speed > 0) {
            double T = (1 - mw.amp - ball.position) / ball.speed;

            auto fmu = [&](double dtnew) {
                return mw.amp * std::cos(dtnew + mw.phi0 + tnow + T)
                    + mw.amp - ball.speed * dtnew;
            };

            
            std::vector<double> roots = find_roots(1e-15, 2 * M_PI, 100, fmu);
            
            
            double tmin = T;

            if (!roots.empty()) {
                // If there are roots, find the smallest positive root
            for (double root : roots) {
                if (root > 0 && root < tmin) {
                     tmin = root;
                    }
                }
                // Add T to tmin if a positive root was found
                tmin += T;
            } 

            return tmin;
        }
        else {
            return std::numeric_limits<double>::infinity();
        }
    }


    double tlwall() {
            double T = -ball.position / ball.speed;
            return (T > 0) ? T : std::numeric_limits<double>::infinity();
        }

    int rwallhit() {
        int col = collision();
        while (col == 2) {
            col = collision();
            if (col == -1) {
                break;
            }
        }
        return col;
    }

    int collision() {
        double trw = trwall();
        double tlw = tlwall();
        // Print debug statements
        // std::cout << "Current time: " << tnow << std::endl;
        // std::cout << "trw, tlw: " << trw << ", " << tlw << std::endl;

        if (tlw < trw) {
            tnow += tlw;
            mw.update(tnow);
            ball.update(-ball.speed, 0);
            
            // Print debug statements
            // std::cout << "Hit left wall" << std::endl;
            // std::cout << "mw.position: " << mw.position << std::endl;
            // std::cout << "ball.position: " << ball.position << std::endl;
            // std::cout << "ball.speed: " << ball.speed << std::endl;
            // std::cout << "mw.speed: " << mw.speed << std::endl;

            return 2;
        } 
        else if (trw < tlw) {
            tnow += trw;
            mw.update(tnow);
            ball.update(2 * mw.speed - ball.speed, mw.position);
            // if (ball.speed < 0) {
            //     ball.update(2 * mw.speed + ball.speed, mw.position);
            // } else {
            //     ball.update(2 * mw.speed - ball.speed, mw.position);
            // }

            //Print debug statements
            // std::cout << "Hit right wall" << std::endl;
            // std::cout << "mw.position: " << mw.position << std::endl;
            // std::cout << "ball.position: " << ball.position << std::endl;
            // std::cout << "ball.speed: " << ball.speed << std::endl;
            // std::cout << "mw.speed: " << mw.speed << std::endl;

            return 1;
        } 
        else {
            std::cout << "Error: Very slow dynamics, similar speed of wall and ball in same direction" << std::endl;
            return -1;
        }
    }

    //calculation via differences
    std::vector<double> difft(int Nv0, int Nph, int Niter) const override {
        std::cout << "WARNING: this calculation is only valid if phase space" 
        "is finite and homogenus (which in this case is not true!)." << '\n';
        std::vector<double> diffarr(100, 0.0);
        auto v0 = linspace(1e-4, 2e-4, Nv0);
        std::vector<double> phases = linspace(0.0, 2 * M_PI, Nph);
        int NP = Nv0 * Nph;

        for (size_t i = 0; i < v0.size(); i++) {
            for (size_t j = 0; j < phases.size(); ++j) {
                FUM2 fermi(ball.mass, -v0[i], mw.amp, phases[j]);
                std::vector<double> varr = {-v0[i]};
                int k = 0;
                for (int iter = 1; iter <= Niter; ++iter) {
                    int hit = fermi.rwallhit();
                    if (hit == -1) {
                        std::cout << -v0[i] << ", " << phases[j] << std::endl;
                        NP -= 1;
                        break;
                    }
                    varr.push_back(fermi.ball.speed);
                    if (iter % 100 == 0 || iter==1) {
                        double D  = diff(varr);
                        double Dcor = diffcor(varr);
                        diffarr[k] += (D)-(1.0/static_cast<double>(iter))*Dcor;
                        k += 1;
                    }
                }
            }
            std::cout << "done: " << static_cast<double>(i + 1) / v0.size() << std::endl;
        }
        std::vector<double> result(100, 0.0);
        for (size_t m = 0; m < diffarr.size(); ++m) {
            result[m] = (1.0 / (NP )) * diffarr[m];
        }
        return result;
    }

    std::vector<double> diffusion(int Nv0, int Nph, int Niter) const override {
        auto v0 = linspace(1e-4, 2e-4, Nv0); 
        auto phases = linspace(0.0, 2 * M_PI, Nph);
        int NP = Nv0 * Nph;
        std::vector<double> narr(100, 1.0);
        double diffu = 0;
        std::vector<double> intermediate_diffusions(100, 0.0);
        double dif = 0;

        for (size_t i = 0; i < v0.size(); i++) {
            for (size_t j = 0; j < phases.size(); j++) {
                FUM2 fermi(1.0, -v0[i],mw.amp, phases[j]); 
                int k = 0;

                for (int iter = 1; iter <= Niter; iter++) {
                    int hit = fermi.rwallhit();
                    if (hit == -1) {
                        std::cout << -v0[i] << ", " << phases[j] << std::endl;
                        NP -= 1;
                        break;
                    }

                    if (iter % 100 == 0 || iter==1) {
                        double speed_diff = std::pow(std::abs(fermi.ball.speed) - (v0[i]), 2);
                        narr[k] = static_cast<double>(iter);
                        intermediate_diffusions[k] += speed_diff;
                        k += 1;
                    }
                }
                diffu += std::pow(std::abs(fermi.ball.speed) - (v0[i]), 2);
            }
            std::cout << "done: " << static_cast<double>(i + 1) / v0.size() << std::endl;
        }
        narr[0]=1;
        for (size_t i = 0; i < intermediate_diffusions.size(); i++)
        {
           intermediate_diffusions[i] /= (narr[i]*2*NP);
        }
        std::cout << (diffu/(Niter*2*NP)) << '\n';
        return intermediate_diffusions;
    }

    
};