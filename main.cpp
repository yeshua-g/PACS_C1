#include "GetPot"
#include <Eigen>

int main(){

    GetPot data("dataGetPot");
    
    if(data("Method/choice","err")==std::string("Gradient")){
        //Gradient descent method
    }
    else if(data("Method/choice","err")==std::string("Heavy-ball")){
       //Heavy_ball method
    }
    else if(data("Method/choice","err")==std::string("Nesterov")){
        //Nesterov method
    }
    else if(data("Method/choice","err")==std::string("Adam")){
        //Adam method
    }else std::cerr << "Declared method not valid!" << std::endl;
 
    return 0;
}