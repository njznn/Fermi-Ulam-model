#include "helper.h"

class GeneralFUM : public Observables {
public:
    //MUST BE CAREFUL ABOUT INITIAL CONDITIONS!!!
    GeneralFUM(const std::vector<double>& marr, const std::vector<double>& x0arr,
     const std::vector<double>& v0arr, double a0, double phi0)
        : tnow(0.0), mw(a0, tnow, phi0), leftwall(0.0) {
        double poss = mw.position;
        if (a0 > std::abs(v0arr.back())){ //&& (phi0 < M_PI/2 || phi0 > 3*M_PI/2)) {
            poss = 1-1.1*a0;
        }
        poss = 1-1.1*a0;
        for (size_t i = 0; i < v0arr.size(); ++i) {
            barr.emplace_back(marr[i], v0arr[i], x0arr[i]);
        }
        barr.back().position = poss;
    }

    double tnow;
    MovingWall mw;
    double leftwall;
    std::vector<Ball> barr;
    double trwall() {
        if (barr.back().position >= (1 - mw.amp)) {
            auto fmu = [&](double dtnew) {
                return mw.amp * std::cos(dtnew + mw.phi0 + tnow)
                    - (barr.back().position - 1) - barr.back().speed * dtnew;
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
        else if (barr.back().speed > 0) {
            double T = (1 - mw.amp - barr.back().position) / barr.back().speed;

            auto fmu = [&](double dtnew) {
                return mw.amp * std::cos(dtnew + mw.phi0 + tnow + T)
                    + mw.amp - barr.back().speed * dtnew;
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
        double T = -barr[0].position / barr[0].speed;
        return (T > 0) ? T : std::numeric_limits<double>::infinity();
    }

    std::vector<double> tbbarr() {
        std::vector<double> tarr;
        for (size_t i = 0; i < barr.size() - 1; ++i) {
            try {
                double T = (barr[i + 1].position - barr[i].position) / (barr[i].speed - barr[i + 1].speed);
                tarr.push_back((T > 0) ? T : std::numeric_limits<double>::infinity());
            } catch (...) {
                tarr.push_back(std::numeric_limits<double>::infinity());
            }
        }
        return tarr;
    }

    int simultcol(double tmin, double t1, double tol){

        if (std::abs(tmin-t1)> tol){
            return 1;
        }
        else{
            return 0; // we have multi collision (at the same time)
                     
        }
    }

    std::vector<int> karr(double tbmin, std::vector<double> &tballs){
        std::vector<int> k; //all balls that hit simultaneously
        for (size_t i = 0; i < tballs.size(); ++i) {
            if (tballs[i] == tbmin) {
                k.push_back(i);
            }
        }
        return k;
    }

    void ballsupd(const std::vector<int>& karray, double t) {
        std::vector<int> combined_karray = karray;

        for (int k : karray) {
            combined_karray.push_back(k + 1);
        }
        for (size_t i = 0; i < barr.size(); ++i) {
            //check if i is in combined_karray
            if (std::find(combined_karray.begin(), 
                combined_karray.end(), i) == combined_karray.end()) {
                barr[i].update(barr[i].speed, barr[i].position + barr[i].speed * t);
            } 
            // Check if i is in karray
            else if (std::find(karray.begin(), karray.end(), i) != karray.end()) {
                if (i + 1 < barr.size()) {
                    ballscol(barr[i], barr[i + 1], t);
                }
            }
        }
    }

    int lwmincol(double &tlw, double & trw, double &tol){
        if (simultcol(tlw, trw, tol)){
            tnow = tlw + tnow;
            mw.update(tnow);
            barr[0].update(-barr[0].speed, 0);
            for (size_t i = 1; i < barr.size(); ++i) {
                barr[i].update(barr[i].speed, barr[i].position + barr[i].speed * tlw);
        }
        }
        else {
            mw.update(tnow+trw);
            tnow = tlw + tnow;
            barr[0].update(-barr[0].speed, 0);
            for (size_t i = 1; i < (barr.size()-1); ++i) {
                barr[i].update(barr[i].speed, barr[i].position + barr[i].speed * tlw); 
            }
            barr.back().update(2 * mw.speed - barr.back().speed, mw.position);
                return -1;
        }
        return -2;
    }

    int tbarrmincol(double &tbmin, double &trw, double &tol, std::vector<double> &tballs){
        std::vector<int> k= karr(tbmin, tballs); //all balls that hit simultaneously
        
        if (simultcol(tbmin, trw, tol)){
            tnow = tnow + tbmin;
            mw.update(tnow);
            ballsupd(k, tbmin);
        }
        else { // between balls hit and right wall hit simultaneously
            mw.update(tnow+trw);
            tnow = tnow + tbmin;
            ballsupd(k, tbmin);
            barr.back().update(2 * mw.speed - barr.back().speed, mw.position);
                return -1;
        }
        return k.back();

    }
    int trwmincol(double &tbmin, double &trw, double &tlw, std::vector<double> &tballs,
            double &tol){
        //just update right wall hit
        tnow += trw;
        mw.update(tnow);
        barr.back().update(2 * mw.speed - barr.back().speed, mw.position);
        for (size_t i = 0; i < barr.size() - 1; ++i) {
            barr[i].update(barr[i].speed, barr[i].position + barr[i].speed * trw);
        }
        return -1;
    }


    int collision() {
    std::vector<double> tballs = tbbarr();
    double tlw = tlwall();
    double trw = trwall();
    double tbmin = *std::min_element(tballs.begin(), tballs.end());
    
    double tol = 1e-7;
    //std::cout << trw << " " << tballs << " " << tlw << std::endl;
    if (std::isinf(std::min(trw, std::min(tbmin, tlw)))){
            return 404;
        }
    else if (tlw < tbmin && tlw < trw) {
        int state = lwmincol(tlw, trw, tol);
        // Uncomment for debugging
        // std::cout << "coll w l wall" << std::endl;
        // std::cout << barr[1].position << std::endl;
        // std::cout << barr[1].speed << std::endl;
        // std::cout << barr[0].position << std::endl;
        // std::cout << barr[0].speed << std::endl;
        // std::cout << mw.position << std::endl;
        // std::cout << mw.speed << std::endl;
        return state;
        }

    else if (tbmin < tlw && tbmin < trw) {
        int state = tbarrmincol(tbmin, trw, tol, tballs);

        // Optionally print values (uncomment if needed)
        // std::cout << "collision between balls" << std::endl;
        // std::cout << barr[1].position << std::endl;
        // std::cout << barr[1].speed << std::endl;
        // std::cout << barr[0].position << std::endl;
        // std::cout << barr[0].speed << std::endl;
        // std::cout << mw.position << std::endl;
        // std::cout << mw.speed << std::endl;
        return state;  // Return vector of indices k
        }

    else if (trw < tlw && trw < *std::min_element(tballs.begin(), tballs.end())) {
        int state= trwmincol(tbmin, trw, tlw, tballs, tol);
    // Uncomment for debugging
    // std::cout << "coll w r wall" << std::endl;
    // std::cout << barr[1].position << std::endl;
    // std::cout << barr[1].speed << std::endl;
    // std::cout << barr[0].position << std::endl;
    // std::cout << barr[0].speed << std::endl;
    // std::cout << mw.position << std::endl;
    // std::cout << mw.speed << std::endl;

    return state;
    }
    // this were all double collisions with rwall and either lwall or multiple between bals colls
    // now all three can also happen: lwcoll, rwcoll, in in-between balls col
    else if (tlw==trw || trw==tbmin || tbmin==tlw || trw==tlw && tlw==tbmin){
        std::cout << "same times!" << '\n';
        if (tlw==trw){
            lwmincol(tlw, trw, tol);
            return -1;
        }
        else if (tlw==tbmin){
            tnow = tnow + tlw;
            ballsupd(karr(tbmin, tballs), tlw);
            mw.update(tnow);
            barr[0].update(-barr[0].speed, 0);
            return -2;
        }
        else if (trw==tbmin){
            int state = tbarrmincol(tbmin, trw, tol, tballs);
            return state;
        }
        else if (trw==tbmin && tbmin==tlw){
            tnow = tnow + tlw;
            ballsupd(karr(tbmin, tballs), tlw);
            mw.update(tnow);
            barr[0].update(-barr[0].speed, 0);
            barr.back().update(2 * mw.speed - barr.back().speed, mw.position);
            return -1;
        }
        else {
            std::cout << "something is wrong at same times" << '\n';
            return (404);
        }
    }
    else {
        std::cout << "error!" << std::endl;
        //exit(1);
        return 404;
        }
    };

    bool rwallhit() {
        int col = collision();
        while (col != -1) {
            col = collision();
            if (col == 404) {
                return false;
            }
        }
        return true;
    }
    bool rballhit() {
        int col = collision();
        while (col != -1 && col !=(barr.size()-2)) {
            col = collision();
            if (col == 404) {
                return false;
            }
        }
        return true;
    }

    std::vector<double> difft(int Nv0, int Nph, int Niter) const override {
        std::cout << "WARNING: this calculation is only valid if phase space" 
        "is finite and homogenus (which in this case is not true!)." << '\n';
        auto v0 = linspace(1e-4, 2e-4, Nv0);
        auto phases = linspace(0.0, 2 * M_PI, Nph);
        int NP = Nv0 * Nph;

        std::vector<double> masses;
        std::vector<double> positions;
        for (const auto& ball : barr) {
            masses.push_back(ball.mass);
            positions.push_back(ball.position);
        }
        std::vector<double> diffarr(100, 0.0);
        for (size_t i = 0; i < v0.size(); i++) {
            for (size_t j = 0; j < phases.size(); ++j) {
                std::vector<double> v0arr(this->barr.size(), -v0[i]);
                GeneralFUM fermi(masses, positions, v0arr,  mw.amp, phases[j]);
                std::vector<double> varr = {-v0[i]};
                int k = 0;
                for (int iter = 1; iter <= Niter; ++iter) {
                    int hit = fermi.rballhit();
                    if (hit == 404) {
                        std::cout << -v0[i] << ", " << phases[j] << std::endl;
                        NP -= 1;
                        break;
                    }
                    varr.push_back(fermi.barr.back().speed);
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
        std::vector<double> masses;
        std::vector<double> positions;
        for (const auto& ball : barr) {
            masses.push_back(ball.mass);
            positions.push_back(ball.position);
        }

        for (size_t i = 0; i < v0.size(); i++) {
            for (size_t j = 0; j < phases.size(); j++) {
                std::vector<double> v0arr(this->barr.size(), -v0[i]);
                GeneralFUM fermi(masses, positions, v0arr,  mw.amp, phases[j]);
                int k = 0;
                double total_diffu = 0;

                for (int iter = 1; iter <= Niter; iter++) {
                    int hit = fermi.rballhit();
                    if (hit == 404) {
                        std::cout << -v0[i] << ", " << phases[j] << std::endl;
                        NP -= 1;
                        break;
                    }

                    if (iter % 100 == 0 || iter==1) {
                        double speed_diff = std::pow(std::abs(fermi.barr.back().speed) -( v0[i]), 2);
                        narr[k] = static_cast<double>(iter);
                        intermediate_diffusions[k] += speed_diff;
                        k += 1;
                    }
                }
                diffu += std::pow((std::abs(fermi.barr.back().speed) - (-v0[i])), 2);
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
