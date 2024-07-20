#include "helper.h"




class DoubleFUM : public Observables {
public:
    DoubleFUM(double massl, double massr, double v0bl, double v0br, double a0, double phi0)
        : tnow(0.0),mw(a0, tnow, phi0), leftwall(0.0), balll(massl, v0bl, 0.20), 
        ballr(massr, v0br, mw.position){
        double poss = mw.position;
        if (a0 > std::abs(v0br)){ //&& (phi0 < M_PI / 2 || phi0 > 3 * M_PI / 2)) {
            poss = 1- 1.1*a0;
        }
        poss = 1- 1.1*a0;
        ballr.update(v0br, poss);
        }

    double tnow; 
    MovingWall mw;
    double leftwall;
    Ball balll;
    Ball ballr;


    double trwall() {
        if (ballr.position >= (1 - mw.amp)) {
            auto fmu = [&](double dtnew) {
                return mw.amp * std::cos(dtnew + mw.phi0 + tnow)
                    - (ballr.position - 1) - ballr.speed * dtnew;
            };

            
            std::vector<double> roots = find_roots(1e-7, 2 * M_PI, 100, fmu);

            double tmin = std::numeric_limits<double>::infinity();
            for (double root : roots) {
                if (root > 0 && root < tmin) {
                    tmin = root;
                    
                }
            }
            
            return tmin;
        }
        // Case 2: ball.speed > 0
        else if (ballr.speed > 0) {
            double T = (1 - mw.amp - ballr.position) / ballr.speed;

            auto fmu = [&](double dtnew) {
                return mw.amp * std::cos(dtnew + mw.phi0 + tnow + T)
                    + mw.amp - ballr.speed * dtnew;
            };

            std::vector<double> roots = find_roots(1e-15, 2 * M_PI, 100, fmu);
            
            
            double tmin = T;
            double minroot = std::numeric_limits<double>::infinity();

            if (!roots.empty()) {
                for (double root : roots) {
                    if (root > 0 && root < minroot) {
                        minroot = root;
                        }
                    }
                tmin = T+minroot;
            } 
            
            return tmin;
        }
        else {
            return std::numeric_limits<double>::infinity();
        }
    }


    double tlwall() {
        double T = -balll.position / balll.speed;
        return (T > 0) ? T : std::numeric_limits<double>::infinity();
    }

    double tbb() {
        try {
            double T = (ballr.position - balll.position) / (balll.speed - ballr.speed);
            return (T > 0 ? T : std::numeric_limits<double>::infinity());
        }
        catch (const std::runtime_error& e) {
            return std::numeric_limits<double>::infinity();
        }
    }
    
    int simultcol(double tmin, double t1, double tol){

        if (std::abs(tmin-t1)> tol){
            return 1;
        }
        else{
            return 0; // we have multi collision (at the same time)
                     
        }
    }

    int collision() {
        double trrw = trwall();
        double tbbc = tbb();
        double tlw = tlwall();

        double tol = 1e-7;
        //std::cout << trrw << " " << tbbc << " " << tlw << std::endl;
        if (std::isinf(std::min(trrw, std::min(tbbc, tlw)))){
            return 404;
        }
        else if (tlw < tbbc && tlw < trrw) {
            /* collision of left ball with left wall */
            if (simultcol(tlw, trrw, tol)){
                tnow = tlw + tnow;
                mw.update(tnow);
                balll.update(-balll.speed, 0);
                ballr.update(ballr.speed, ballr.position + ballr.speed * tlw);
            }
            else{
                mw.update(tnow+trrw);
                tnow = tlw+tnow;
                balll.update(-balll.speed, 0);
                ballr.update(2 * mw.speed - ballr.speed, mw.position);
                // if (ballr.speed < 0) {
                //     ballr.update(2 * mw.speed + ballr.speed, mw.position);
                // } else {
                //     ballr.update(2 * mw.speed - ballr.speed, mw.position);
                // }
                return -1;
            }
            // std::cout << "coll w l wall" << std::endl;
            // std::cout << ballr.speed << " " << ballr.position << std::endl;
            // std::cout << mw.position << " " << mw.speed << std::endl;
            // std::cout << balll.position << " " << balll.speed << std::endl;
            return 0;
        }
        else if (tbbc < tlw && tbbc < trrw) {
            /* collision between balls */
            if (simultcol(tbbc,trrw, tol)){
                tnow = tbbc + tnow;
                mw.update(tnow);
                ballscol(balll, ballr, tbbc);
                
            }
            else {
                mw.update(tnow+trrw);
                tnow = tbbc + tnow;
                ballscol(balll, ballr, tbbc);
                ballr.update(2 * mw.speed - ballr.speed, mw.position);
                return -1;
            }
           
            // std::cout << "collision between balls" << std::endl;
            // std::cout << ballr.position << std::endl;
            // std::cout << balll.position << std::endl;
            // std::cout << mw.speed << '\n';
            // std::cout << ballr.speed << std::endl;
            // std::cout << balll.speed << std::endl;

            return 1;
        }
        else if (trrw < tlw && trrw < tbbc) {
            /* collision of right ball with right wall */ 
            //first if not possible in this case
            if (simultcol(trrw, tlw, tol) && simultcol(trrw, tbbc, tol)){
                tnow = trrw + tnow;
                mw.update(tnow);
                ballr.update(2 * mw.speed - ballr.speed, mw.position);
                balll.update(balll.speed, balll.position + balll.speed * trrw);
            }
            else if (simultcol(trrw, tlw, tol)==0){
                tnow = trrw+tnow;
                mw.update(tnow);
                balll.update(-balll.speed, 0);
                ballr.update(2 * mw.speed - ballr.speed, mw.position);
                
            }
            else {
                tnow = trrw + tnow;
                ballscol(balll, ballr, tbbc);
                mw.update(tnow);
                ballr.update(2 * mw.speed - ballr.speed, mw.position);
             }

            // std::cout << "rwcoll" << '\n';
            // std::cout << mw.position << std::endl;
            // std::cout << ballr.position << std::endl;
            // std::cout << balll.position << std::endl;
            // std::cout << mw.speed << std::endl;
            // std::cout << ballr.speed << std::endl;
            // std::cout << balll.speed << std::endl;
            return -1;
        }
        else if (tlw== tbbc){
            ballscol(balll, ballr, tbbc);
            balll.update(-balll.speed, 0);
            std::cout << "bb and lw coll simult" << '\n';
            return 0;
        }
        else {
            std::cout << "error, most probably, initial condition problem (wall faster than ball)" << std::endl;
            return 404;
        }
    }
    bool rwallhit(){
        int col = collision();
            while (col != -1 ) {
                col = collision();
                if (col == 404) {
                    return false; 
                }
            }
            return true;
    }
    bool rballhit() {
        int col = collision();
        while (col ==0 ) {
            col = collision();
            if (col == 404) {
                return false; 
            }
        }
        //std::cout << col << '\n';
        return true;
    }
    

    std::vector<double> difft(int Nv0, int Nph, int Niter) const override {
        std::cout << "WARNING: this calculation is only valid if phase space" 
        "is finite and homogenus (which in this case is not true!)." << '\n';
        std::vector<double> diffarr(100, 0.0);
        auto v0 = linspace(1e-4, 2e-4, Nv0);
        std::vector<double> phases = linspace(0.0, 2 * M_PI, Nph);
        int NP = Nv0 * Nph;
        for (size_t i = 0; i < v0.size(); i++) {
            for (size_t j = 0; j < phases.size(); ++j) {
                DoubleFUM fermi(balll.mass,ballr.mass, -v0[i],-v0[i],  mw.amp, phases[j]);
                std::vector<double> varr = {-v0[i]};
                int k = 0;
                for (int iter = 1; iter <= Niter; ++iter) {
                    int hit = fermi.rballhit();
                    if (hit == 404) {
                        std::cout << -v0[i] << ", " << phases[j] << std::endl;
                        NP -= 1;
                        break;
                    }
                    varr.push_back(fermi.ballr.speed);
                    if (iter % 100 == 0 || iter==1) {
                        double D  = diff(varr);
                        double Dcor = diffcor(varr);
                        diffarr[k] += (D-(1.0/static_cast<double>(iter))*Dcor);
                        k += 1;
                    }
                }
            }
            std::cout << "done: " << static_cast<double>(i + 1) / v0.size() << std::endl;
        }
        
        std::vector<double> result(100, 0.0);
        for (size_t m = 0; m < diffarr.size(); ++m) {
            result[m] = (1.0 / (NP)) * diffarr[m];
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
        

        for (size_t i = 0; i < v0.size(); i++) {
            for (size_t j = 0; j < phases.size(); j++) {
                DoubleFUM fermi(balll.mass,ballr.mass, -v0[i],-v0[i],  mw.amp, phases[j]);
                int k = 0;
                double total_diffu = 0;
                double del = 0;

                for (int iter = 1; iter <= Niter; iter++) {
                    int hit = fermi.rballhit();
                    if (hit == 404) {
                        std::cout << -v0[i] << ", " << phases[j] << std::endl;
                        NP -= 1;
                        break;
                    }
                    if (iter % 100 == 0 || iter==1) {
                        double speed_diff = std::pow(std::abs(fermi.ballr.speed) -(v0[i]), 2);
                        narr[k] = static_cast<double>(iter);
                        intermediate_diffusions[k] += speed_diff;
                        k += 1;
                    }
                }
                diffu += std::pow((std::abs(fermi.ballr.speed) - (v0[i])), 2);
            }
            std::cout << "done: " << static_cast<double>(i + 1) / v0.size() << std::endl;
        }
        narr[0]=1;
        
        for (size_t i = 0; i < intermediate_diffusions.size(); i++)
        {
           intermediate_diffusions[i] /= (narr[i]*2*NP);
        }
        return intermediate_diffusions;
    }

};
