#include "Option.hpp"
#include <algorithm>
#include <iostream>
class BasketOption : public Option
{
public:
    double strike_; 
    PnlVect* weights_; //vect poids
    
    BasketOption(double T, int nbTimeSteps, int size, double strike, const PnlVect* weights)
        : Option(T, nbTimeSteps, size), strike_(strike)
    {
        weights_ = pnl_vect_copy(weights); 
    }

    //destructeur    
    ~BasketOption()
    {
        pnl_vect_free(&weights_);
    }

    //payoff
    double payoff(const PnlMat* path) override
    {
        double weightedSum = 0.0;
        for (int d = 0; d < size_; ++d)
        {
            weightedSum += pnl_vect_get(weights_, d) * pnl_mat_get(path, path->m - 1, d);
        }

        double payoff = weightedSum - strike_;
        return std::max(0.0, payoff); 
    }
};