#include "json_helper.hpp"
#include "BlackScholesModel.hpp"
#include "Option.hpp"
#include "AsianOption.cpp"
#include "BasketOption.cpp"
#include "PerformanceOption.cpp"
#include "MonteCarlo.hpp"
#include "PricingResults.hpp"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Wrong number of arguments. Exactly one argument is required" << std::endl;
        std::exit(0);
    }
    std::ifstream ifs(argv[1]);
    nlohmann::json j = nlohmann::json::parse(ifs);
    int size;
    std::string option_type;

    PnlVect* volatility = pnl_vect_new();
    PnlVect* weights = pnl_vect_new();
    PnlVect* spots = pnl_vect_new();

    option_type = j.at("option type").get<std::string>();

    j.at("option size").get_to(size);

    j.at("volatility").get_to(volatility);
    if (volatility->size == 1 && size > 1) {
        pnl_vect_resize_from_scalar(volatility, size, GET(volatility, 0));
    }
    
    j.at("spot").get_to(spots);
    if (spots->size == 1 && size > 1) {
        pnl_vect_resize_from_scalar(spots, size, GET(spots, 0));
    }

    j.at("payoff coefficients").get_to(weights);
    if (weights->size == 1 && size > 1) {
        pnl_vect_resize_from_scalar(weights, size, GET(weights, 0));
    }
    
    double maturity;
    j.at("maturity").get_to(maturity);

    double strike;

    double r;
    j.at("interest rate").get_to(r);

    double correlation;
    j.at("correlation").get_to(correlation);

    int nbTimeSteps;
    j.at("timestep number").get_to(nbTimeSteps);

    int nbSamples;
    j.at("sample number").get_to(nbSamples);

    double fdStep;
    j.at("fd step").get_to(fdStep);

    PnlRng* rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    BlackScholesModel* blackScholesModel = new BlackScholesModel(size, r, correlation, volatility, spots);;

    Option* option = nullptr;

    if (option_type == "basket") {
        j.at("strike").get_to(strike);
        option = new BasketOption(maturity, nbTimeSteps, size, strike, weights);
    } else if (option_type == "performance") {
        option = new PerformanceOption(maturity, nbTimeSteps, size, weights);
    } else if (option_type == "asian") {
        j.at("strike").get_to(strike);
        option = new AsianOption(maturity, nbTimeSteps, size, strike, weights);
    }

    MonteCarlo* monteCarlo = new MonteCarlo(blackScholesModel, option, rng, fdStep, nbSamples);

    double price = 0.0;
    double priceStdDev = 0.0;
    PnlVect* delta = pnl_vect_create(size);
    PnlVect* deltaStdDev = pnl_vect_create(size);

    monteCarlo->deltaPrice(price, priceStdDev, delta, deltaStdDev);

    PricingResults res(price, priceStdDev, delta, deltaStdDev);
    std::cout << res << std::endl;

    pnl_vect_free(&delta);
    pnl_vect_free(&deltaStdDev);
    pnl_vect_free(&spots);
    pnl_vect_free(&weights);
    pnl_vect_free(&volatility);
    pnl_rng_free(&rng);
    exit(0);
}