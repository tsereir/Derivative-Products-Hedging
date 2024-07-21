#include"Option.hpp"
class PerformanceOption : public Option
{
    public:
    PnlVect* weights_; //vect poids
    PerformanceOption(double T, int nbTimeSteps, int size, const PnlVect* weights)
        : Option(T, nbTimeSteps, size), weights_(pnl_vect_copy(weights))
    {
    }
    ~PerformanceOption()
    {
        pnl_vect_free(&weights_);
    }
    double payoff(const PnlMat* path) override
    {

        double performance = 1.0;
        for (int i=1; i< nbTimeSteps_ + 1; i++)
        {
            //somme des prix a l'intant ti
            double sum_ti = 0.0;
            double sum_ti_minus_1 = 0.0; //Somme  à l'instant ti-1
        
            for(int d= 0; d < size_; d++)
            {
                //somme pondérees des 2 instants
                sum_ti += pnl_mat_get(path, i, d) * pnl_vect_get(weights_, d);
                sum_ti_minus_1 += pnl_mat_get(path, i-1, d) * pnl_vect_get(weights_, d);
            }
            //terme de performance a l'intant ti, c'ets le rendement relatif entre les 2 instants
            double term = (sum_ti / sum_ti_minus_1) - 1.0;
            term = (term > 0.0) ? term : 0.0;
            //maj de la performance 
            performance += term;   
        }
        return performance ;   
    }
};