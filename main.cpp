#include "minimization.hpp"

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

    else if(data.Method=="HeavyBall"){
    switch (data.step_size_strategy) {
        case StepSizeStrategy::ExponentialDecay:
            res = heavy_ball<StepSizeStrategy::ExponentialDecay>(data);
            break;
        case StepSizeStrategy::InverseDecay:
            res = heavy_ball<StepSizeStrategy::InverseDecay>(data);
            break;
         case StepSizeStrategy::FixedAlpha:
            res = heavy_ball<StepSizeStrategy::FixedAlpha>(data);
            break;
    }
    }

    else if(data.Method=="Nesterov"){
    switch (data.step_size_strategy) {
        case StepSizeStrategy::ExponentialDecay:
            res = Nesterov<StepSizeStrategy::ExponentialDecay>(data);
            break;
        case StepSizeStrategy::InverseDecay:
            res = Nesterov<StepSizeStrategy::InverseDecay>(data);
            break;
        case StepSizeStrategy::FixedAlpha:
            res = Nesterov<StepSizeStrategy::FixedAlpha>(data);
            break;
    }
    }

    else if(data.Method=="Adam"){
        res=Adam(data);
    }

    else throw std::invalid_argument("Invalid method");



    print(res);
    

    return 0;
}

