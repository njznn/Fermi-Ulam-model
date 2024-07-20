
#include "helper.h"




std::vector<double> linspace(double start, double end, int num) {
    if (num <= 0) {
        throw std::invalid_argument("Number of elements must be greater than zero.");
    } else if (num == 1) {
        return {start}; // Return a vector with just the start value
    }

    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        result[i] = start + step * i;
    }
    return result;
}


double diff(const std::vector<double>& arr) {
        double res = 0;
        double w0 = std::abs(arr[1]) -std::abs(arr[0]);
        for (size_t i = 0; i < arr.size() - 1; ++i) {
            if (i==0){
                res += (std::abs(arr[i + 1]) - std::abs(arr[i]))*w0*0.5;
            }
            else {
                res += (std::abs(arr[i + 1]) - std::abs(arr[i]))*w0;
            }

        }
        return res;
    }
    
double diffcor(const std::vector<double>& arr) {
        double res = 0;
        double w0 =  std::abs(arr[1]) - std::abs(arr[0]);
        for (size_t i = 0; i < arr.size() - 1; ++i) {
            if (i==0){
                res += (std::abs(arr[i + 1]) - std::abs(arr[i]))*w0*static_cast<double>(i)*0.5;
            }
            else {
                res += (std::abs(arr[i + 1]) - std::abs(arr[i]))*w0*static_cast<double>(i);;
            }

        }
        return res;
    }


    Ball::Ball(double mass, double speed, double position)
        : mass(mass), speed(speed), position(position) {}

    void Ball::update(double newSpeed, double newPosition) {
        speed = newSpeed;
        position = newPosition;
    }

// CORRECT WALL POSIITON

    MovingWall::MovingWall(double amp, double t, double phi0)
        : amp(amp), phi0(phi0) {
        update(t);
    }

    
    void MovingWall::update(double tnew) {
        position = 1+amp * std::cos(tnew + phi0);
        speed = -amp * std::sin(tnew + phi0);
        phi = std::fmod(tnew + phi0, 2 * M_PI);
    }


void ballscol(Ball& ball1, Ball& ball2, double tcol) {
    // Calculate new velocities
    double v1 = (((ball1.mass - ball2.mass) / (ball1.mass + ball2.mass)) * ball1.speed +
                 (2 * ball2.mass / (ball1.mass + ball2.mass)) * ball2.speed);

    double v2 = (((ball2.mass - ball1.mass) / (ball1.mass + ball2.mass)) * ball2.speed +
                 (2 * ball1.mass / (ball1.mass + ball2.mass)) * ball1.speed);

    
    double xcol = ball1.position + ball1.speed * tcol;

    // Update the balls' states
    ball1.update(v1, xcol);
    ball2.update(v2, xcol);
}

std::vector<double> find_roots(double xmin, double xmax, int N, std::function<double(double)> f) {
    std::vector<double> roots;

    // Generate evenly spaced x values
    std::vector<double> x_range;
    double step = (xmax - xmin) / (N - 1);
    for (int i = 0; i < N; ++i) {
        x_range.push_back(xmin + i * step);
    }
    boost::uintmax_t max_iter = 1e10;
    const int significant_digits = 10;
    // Iterate over the range and find roots
    for (size_t i = 0; i < x_range.size() - 1; ++i) {
        double x1 = x_range[i];
        double x2 = x_range[i + 1];
        double f1 = f(x1);
        double f2 = f(x2);

        // Check for sign change
        if (f1 * f2 < 0) {
            // Use Boost's brent_find_minima to find root within [x1, x2]
            
            auto result = boost::math::tools::toms748_solve(f, x1, x2,
             boost::math::tools::eps_tolerance<double>(significant_digits),max_iter);
            std::streamsize precision_1 = std::cout.precision(std::numeric_limits<double>::digits10);
            if (result.first) {
                roots.push_back(result.first);
            }
        }
    }

    return roots;
}

Observables::~Observables() {}