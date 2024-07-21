#include "MonteCarlo.hpp"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include <cmath>
#include <iostream>
MonteCarlo::MonteCarlo(BlackScholesModel* mod, Option* opt, PnlRng* rng, double fdStep, long nbSamples)
        : mod_(mod), opt_(opt), rng_(rng), fdStep_(fdStep), nbSamples_(nbSamples) 
{
}


void MonteCarlo::deltaPrice(double& prix, double& std, PnlVect* delta, PnlVect* std_dev){
    PnlMat* path = pnl_mat_create(opt_->nbTimeSteps_ + 1, opt_->size_);
    double var_price = 0;
    double dfPayOff = 0;
    for (size_t i = 0; i < nbSamples_; i++){
        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
        double payOff = opt_->payoff(path);
        prix += payOff;
        var_price += payOff * payOff;
        for (size_t d = 0; d < opt_->size_; d++){
            mod_->shiftAsset(path, path, d, fdStep_, 0, opt_->nbTimeSteps_);
            dfPayOff = opt_->payoff(path);
            mod_->shiftAsset(path, path, d, (-2 * fdStep_) / (1 + fdStep_), 0, opt_->nbTimeSteps_);
            dfPayOff -= opt_->payoff(path);
            mod_->shiftAsset(path, path, d, fdStep_ / (1 - fdStep_), 0, opt_->nbTimeSteps_);
            dfPayOff /= (2 * fdStep_ * MGET(path, 0, d));
            LET(delta, d) = GET(delta, d) +  dfPayOff;
            LET(std_dev, d) = GET(std_dev, d) +  dfPayOff * dfPayOff;
            dfPayOff = 0;
        }
    }
    prix = (prix / nbSamples_);
    var_price = exp(-2 * mod_->r_ * opt_->T_) * ((var_price / nbSamples_) - prix * prix);
    prix *= exp(-mod_->r_ * opt_->T_);
    std = sqrt(var_price/nbSamples_);
    for (size_t d = 0; d < opt_->size_; d++){
        LET(delta, d) = exp(-mod_->r_ * opt_->T_) * (GET(delta, d) / nbSamples_);
        LET(std_dev, d) = exp(-2 * mod_->r_ * opt_->T_) * ((GET(std_dev, d) / nbSamples_) - GET(delta, d) * GET(delta, d));
        LET(std_dev, d) = sqrt(GET(std_dev, d)/nbSamples_);
    }
    pnl_mat_free(&path);
}

void MonteCarlo::deltaPrice(PnlMat* path, const PnlMat* past, double t, double& prix, double& std, PnlVect* delta, PnlVect* std_dev){
    double timeStep = opt_->T_ / opt_->nbTimeSteps_;
    size_t iPlus1 = past->m - 1;
    PnlVect* s_t = pnl_vect_create(opt_->size_);
    pnl_mat_get_row(s_t,past,iPlus1);
    double var_price = 0;
    double dfPayOff = 0;
    for (size_t i = 0; i < nbSamples_; i++){
        mod_->asset(path, t, opt_->T_, opt_->nbTimeSteps_, rng_, past);
        for (size_t d = 0; d < opt_->size_; d++){
            mod_->shiftAsset(path, path, d, fdStep_, t, opt_->nbTimeSteps_);
            dfPayOff = opt_->payoff(path);
            mod_->shiftAsset(path, path, d, (-2 * fdStep_) / (1 + fdStep_), t, opt_->nbTimeSteps_);
            dfPayOff -= opt_->payoff(path);
            mod_->shiftAsset(path, path, d, fdStep_ / (1 - fdStep_), t, opt_->nbTimeSteps_);
            LET(delta, d) = GET(delta, d) +  dfPayOff;
            dfPayOff = 0;
        }
    }
    for (size_t d = 0; d < opt_->size_; d++){
        LET(delta, d) = exp(-mod_->r_ * (opt_->T_ - t)) * (GET(delta, d) / (nbSamples_ * (2 * fdStep_ * GET(s_t, d))));
    }
    pnl_vect_free(&s_t);
}
