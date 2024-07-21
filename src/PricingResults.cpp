#include <iostream>
#include "PricingResults.hpp"
#include "json_helper.hpp"

std::ostream&
operator<<(std::ostream& stm, const PricingResults& res)
{

    nlohmann::json j = {
        {"price", res.price},
        {"priceStdDev", res.priceStdDev},
        {"delta", res.delta},
        {"deltaStdDev", res.deltaStdDev}
    };
    stm << j << std::endl;
    return stm;
}