#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include "GA/c3ga.h"
#include <vector>

c3ga::rotor GARotorEstimator(const std::vector<Eigen::Vector3d>& P, const std::vector<Eigen::Vector3d>& Q, const std::vector<double>& w);
c3ga::rotor GARotorEstimator(const Eigen::Matrix3d& M, const double S);
