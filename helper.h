#ifndef HELPER_H
#define HELPER_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <fstream>
//#include <Eigen/Dense>
#include <boost/math/tools/roots.hpp>



constexpr double INF = std::numeric_limits<double>::infinity();

std::vector<double> linspace(double start, double end, int num);
double diff(const std::vector<double>& arr);
double diffcor(const std::vector<double>& arr);

class Ball {
public:
    Ball(double mass, double speed, double position);
    void update(double newSpeed, double newPosition);
    double mass;
    double speed;
    double position;
};


class MovingWall {
public:
    MovingWall(double amp, double t, double phi0);
    void update(double tnew);

    double amp;
    double phi0;
    double position;
    double speed;
    double phi;
};

void ballscol(Ball& ball1, Ball& ball2, double tcol);

std::vector<double> find_roots(double xmin, double xmax, int N, std::function<double(double)> f);

class Observables {
public:
    virtual ~Observables();

    virtual std::vector<double> diffusion(int Nv0, int Nph, int Niter) const=0;
    virtual std::vector<double> difft(int Nv0, int Nph, int Niter) const=0;
};
#endif // HELPER_H