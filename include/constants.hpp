#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <vector>

namespace vkl {

  class Constants{
  public:
    // Finite difference coefficients and relative indices.
    // 1st derivatives
    static const std::vector<int>    derivative_1_forward_1_index;
    static const std::vector<double> derivative_1_forward_1_coeff;
    static const std::vector<int>    derivative_1_forward_2_index;
    static const std::vector<double> derivative_1_forward_2_coeff;
    static const std::vector<int>    derivative_1_backward_1_index;
    static const std::vector<double> derivative_1_backward_1_coeff;
    static const std::vector<int>    derivative_1_backward_2_index;
    static const std::vector<double> derivative_1_backward_2_coeff;
    static const std::vector<int>    derivative_1_central_2_index;
    static const std::vector<double> derivative_1_central_2_coeff;
    // 2nd derivatives
    static const std::vector<int>    derivative_2_forward_1_index;
    static const std::vector<double> derivative_2_forward_1_coeff;
    static const std::vector<int>    derivative_2_forward_2_index;
    static const std::vector<double> derivative_2_forward_2_coeff;
    static const std::vector<int>    derivative_2_backward_1_index;
    static const std::vector<double> derivative_2_backward_1_coeff;
    static const std::vector<int>    derivative_2_backward_2_index;
    static const std::vector<double> derivative_2_backward_2_coeff;
    static const std::vector<int>    derivative_2_central_2_index;
    static const std::vector<double> derivative_2_central_2_coeff;
  };

}

#endif /* CONSTANTS_HPP */
