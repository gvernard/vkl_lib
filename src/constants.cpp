#include "constants.hpp"

using namespace vkl;


// Definition of static variables:
// 2nd derivatives
const std::vector<int>    Constants::derivative_1_forward_1_index({0,1});
const std::vector<double> Constants::derivative_1_forward_1_coeff({-1,1});
const std::vector<int>    Constants::derivative_1_forward_2_index({0,1,2});
const std::vector<double> Constants::derivative_1_forward_2_coeff({-1.5,2.0,-0.5});
const std::vector<int>    Constants::derivative_1_backward_1_index({-1,0});
const std::vector<double> Constants::derivative_1_backward_1_coeff({-1,1});
const std::vector<int>    Constants::derivative_1_backward_2_index({-2,-1,0});
const std::vector<double> Constants::derivative_1_backward_2_coeff({0.5,-2,1.5});
const std::vector<int>    Constants::derivative_1_central_2_index({-1,0,1});
const std::vector<double> Constants::derivative_1_central_2_coeff({-0.5,0.0,0.5});
// 2nd derivatives
const std::vector<int>    Constants::derivative_2_forward_1_index({0,1,2});
const std::vector<double> Constants::derivative_2_forward_1_coeff({1,-2,1});  
const std::vector<int>    Constants::derivative_2_forward_2_index({0,1,2,3});
const std::vector<double> Constants::derivative_2_forward_2_coeff({2,-5,4,-1});
const std::vector<int>    Constants::derivative_2_backward_1_index({-2,-1,0});
const std::vector<double> Constants::derivative_2_backward_1_coeff({1,-2,1});  
const std::vector<int>    Constants::derivative_2_backward_2_index({-3,-2,-1,0});
const std::vector<double> Constants::derivative_2_backward_2_coeff({-1,4,-5,2});
const std::vector<int>    Constants::derivative_2_central_2_index({-1,0,1});
const std::vector<double> Constants::derivative_2_central_2_coeff({1,-2,1});
