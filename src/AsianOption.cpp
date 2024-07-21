#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

class AsianOption : public Option
{
private:
    double strike_;               // Prix d'exercice de l'option 
    PnlVect* weights_; // Poids des sous-jacents dans le calcul de la moyenne

public:
    AsianOption(double T, int nbTimeSteps, int size, double strike, const PnlVect* weights)
        :Option(T, nbTimeSteps, size), strike_(strike), weights_(pnl_vect_copy(weights))
    {
        
    }

    double payoff(const PnlMat* path) override
    { 
        PnlVect* average_prices = pnl_vect_create_from_scalar(size_, 0.0);
        for (int j = 0; j < size_; j++)
        {
            double weighted_sum = 0.0;
            double x = 0.0;
            for (int i = 0; i < nbTimeSteps_+1; i++)
            {
                x = MGET(path,i,j);
                weighted_sum += MGET(path, i, j);
            }
            LET(average_prices, j) = (GET(weights_, j) *weighted_sum) / (nbTimeSteps_ + 1);
        }
        double option_payoff = std::max(pnl_vect_sum(average_prices) - strike_, 0.0);
        pnl_vect_free(&average_prices);
        return option_payoff;
    }
};
