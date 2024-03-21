#include "gradient.hpp"

int main(){
    OptimizationParameters data=read_optimization_parameters("dataGetPot");

    std::vector<double> res(data.dimension);

    switch (data.step_size_strategy) {
        case StepSizeStrategy::ExponentialDecay:
            res = gradient_descent<StepSizeStrategy::ExponentialDecay>(data);
            break;
        case StepSizeStrategy::InverseDecay:
            res = gradient_descent<StepSizeStrategy::InverseDecay>(data);
            break;
        case StepSizeStrategy::ApproximateLineSearch:
            res = gradient_descent<StepSizeStrategy::ApproximateLineSearch>(data);
            break;
    }


    print(res);
    

    return 0;
}

