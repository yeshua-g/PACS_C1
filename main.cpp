#include "gradient.hpp"

int main(){
    OptimizationParameters data=read_optimization_parameters("dataGetPot");

    std::vector<double> res(data.dimension);

    if(data.Method=="Gradient"){
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
    }

    if(data.Method=="HeavyBall"){
    switch (data.step_size_strategy) {
        case StepSizeStrategy::ExponentialDecay:
            res = heavy_ball<StepSizeStrategy::ExponentialDecay>(data);
            break;
        case StepSizeStrategy::InverseDecay:
            res = heavy_ball<StepSizeStrategy::InverseDecay>(data);
            break;
    }
    }


    print(res);
    

    return 0;
}

