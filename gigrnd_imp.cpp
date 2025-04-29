#include <cmath>
#include <random>
#include <pybind11/pybind11.h>

namespace py = pybind11;

// static random number generator at file scope
static std::mt19937 gen(42);

// set the seed
void set_seed(unsigned int seed) {
    gen.seed(seed);
}

// psi
double psi(double x, double alpha, double lam) {
    return -alpha * (std::cosh(x) - 1.0) - lam * (std::exp(x) - x - 1.0);
}

// dpsi
double dpsi(double x, double alpha, double lam) {
    return -alpha * std::sinh(x) - lam * (std::exp(x) - 1.0);
}

// g
double g(double x, double sd, double td, double f1, double f2) {
    if (x >= -sd && x <= td) {
        return 1.0;
    }
    else if (x > td) {
        return f1;
    }
    else { // x < -sd
        return f2;
    }
}

// Main function to generate GIG random variate
double gigrnd(double p, double a, double b) {
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Initial parameter setup
    double lam = p;
    double omega = std::sqrt(a * b);

    // Handle negative lam
    bool swap = false;
    if (lam < 0) {
        lam = -lam;
        swap = true;
    }

    // Compute alpha
    double alpha = std::sqrt(std::pow(omega, 2) + std::pow(lam, 2)) - lam;

    // Find t
    double x = -psi(1.0, alpha, lam);
    double t;
    if (x >= 0.5 && x <= 2.0) {
        t = 1.0;
    }
    else if (x > 2.0) {
        if (alpha == 0.0 && lam == 0.0) {
            t = 1.0;
        }
        else {
            t = std::sqrt(2.0 / (alpha + lam));
        }
    }
    else {  // x < 0.5
        if (alpha == 0.0 && lam == 0.0) {
            t = 1.0;
        }
        else {
            t = std::log(4.0 / (alpha + 2.0 * lam));
        }
    }

    // Find s
    x = -psi(-1.0, alpha, lam);
    double s;
    if (x >= 0.5 && x <= 2.0) {
        s = 1.0;
    }
    else if (x > 2.0) {
        if (alpha == 0.0 && lam == 0.0) {
            s = 1.0;
        }
        else {
            s = std::sqrt(4.0 / (alpha * std::cosh(1.0) + lam));
        }
    }
    else {  // x < 0.5
        if (alpha == 0.0 && lam == 0.0) {
            s = 1.0;
        }
        else if (alpha == 0.0) {
            s = 1.0 / lam;
        }
        else if (lam == 0.0) {
            s = std::log(1.0 + 1.0 / alpha + std::sqrt(1.0 / std::pow(alpha, 2) + 2.0 / alpha));
        }
        else {
            double temp1 = 1.0 / lam;
            double temp2 = std::log(1.0 + 1.0 / alpha + std::sqrt(1.0 / std::pow(alpha, 2) + 2.0 / alpha));
            s = std::min(temp1, temp2);
        }
    }

    double eta = -psi(t, alpha, lam);
    double zeta = -dpsi(t, alpha, lam);
    double theta = -psi(-s, alpha, lam);
    double xi = dpsi(-s, alpha, lam);

    double p_aux = 1.0 / xi;
    double r = 1.0 / zeta;

    double td = t - r * eta;
    double sd = s - p_aux * theta;
    double q = td + sd;

    // main sampling loop
    double rnd;
    while (true) {
        double U = dis(gen);
        double V = dis(gen);
        double W = dis(gen);
        if (U < q / (p_aux + q + r)) {
            rnd = -sd + q * V;
        }
        else if (U < (q + r) / (p_aux + q + r)) {
            rnd = td - r * std::log(V);
        }
        else {
            rnd = -sd + p_aux * std::log(V);
        }

        double f1 = std::exp(-eta - zeta * (rnd - t));
        double f2 = std::exp(-theta + xi * (rnd + s));
        double g_val = g(rnd, sd, td, f1, f2);
        double psi_val = psi(rnd, alpha, lam);
        if (W * g_val <= std::exp(psi_val)) {
            break;
        }
    }

    rnd = std::exp(rnd) * (lam / omega + std::sqrt(1.0 + std::pow(lam, 2) / std::pow(omega, 2)));
    if (swap) {
        rnd = 1.0 / rnd;
    }
    rnd = rnd / std::sqrt(a / b);

    return rnd;
}

PYBIND11_MODULE(gigrnd_imp, m) {
    m.doc() = "C++ implementation of Generalized Inverse Gaussian random number generator";

    m.def("set_seed", &set_seed, "Set the random number generator seed");
    m.def("psi", &psi, "GIG psi function");
    m.def("dpsi", &dpsi, "GIG dpsi function");
    m.def("g", &g, "GIG g function");
    m.def("gigrnd", &gigrnd, "Generate a random variate from the Generalized Inverse Gaussian distribution",
        py::arg("p"), py::arg("a"), py::arg("b"));

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}