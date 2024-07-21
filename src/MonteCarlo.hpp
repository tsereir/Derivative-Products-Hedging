#pragma once

#include "Option.hpp"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_random.h"

class MonteCarlo
{
  public:
    BlackScholesModel* mod_; /*! pointeur vers le modèle */
    Option* opt_;            /*! pointeur sur l'option */
    PnlRng* rng_;            /*! pointeur sur le générateur */
    double fdStep_;          /*! pas de différence finie */
    long nbSamples_;         /*! nombre de tirages Monte Carlo */


    /**
     * Declaration de Constructeur
     */
    MonteCarlo(BlackScholesModel* mod, Option* opt, PnlRng* rng, double fdStep, long nbSamples);

    void deltaPrice(double& prix, double& std, PnlVect* delta, PnlVect* std_dev);
    
    void deltaPrice(PnlMat* path, const PnlMat* past, double t, double& prix, double& std, PnlVect* delta, PnlVect* std_dev);
};
