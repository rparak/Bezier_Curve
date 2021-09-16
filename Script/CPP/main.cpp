/**
 * MIT License
 * Copyright(c) 2021 Roman Parak
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * Author   : Roman Parak
 * Email    : Roman.Parak@outlook.com
 * File Name: main.cpp
 * Github   : https://github.com/rparak
 */

// Own library for Bézier curve calculation
#include "Bezier.h"

int main(){
    /*
     Note:
        $ g++ -o test -std=c++17 main.cpp
        $ ./test
    */

    // Set input parameters
    /*
    std::vector<std::vector<float>> points_lin = {{1.0, 1.0, 1.0}, {1.25, 2.0, 2.5}};
    std::vector<std::vector<float>> points_quadratic = {{1.0, 1.0, 1.0}, {1.25, 2.0, 2.5}, 
                                                        {1.75, 2.0, 1.5}};
    std::vector<std::vector<float>> points_cubic = {{1.0, 1.0, 1.0}, {1.25, 2.0, 2.5}, 
                                                    {1.75, 2.0, 1.5}, {2.0, 1.0, 1.0}};
    std::vector<std::vector<float>> points_n_deg = {{1.0, 1.0, 1.0}, {1.25, 2.0, 2.5}, 
                                                    {1.75, 2.0, 1.5}, {2.0, 1.0, 1.0}, 
                                                    {1.0, -1.0, 2.0}, {1.25, -2.0, 1.75}, 
                                                    {1.75, -2.0, 2.75}, {2.0, -1.0, 2.0}};
    */

    std::vector<std::vector<float>> points_lin = {{1.0, 1.0}, {1.25, 2.0}};
    std::vector<std::vector<float>> points_quadratic = {{1.0, 1.0}, {1.25, 2.0}, 
                                                        {1.75, 2.0}};
    std::vector<std::vector<float>> points_cubic = {{1.0, 1.0}, {1.25, 2.0}, 
                                                    {1.75, 2.0}, {2.0, 1.0}};
    std::vector<std::vector<float>> points_n_deg = {{1.0, 1.0}, {1.25, 2.0}, 
                                                    {1.75, 2.0}, {2.0, 1.0}, 
                                                    {1.0, -1.0}, {1.25, -2.0}, 
                                                    {1.75, -2.0}, {2.0, -1.0}};

    // Display curve points to the console
    bool display_curve = false;

    // Number of samples to generate. Must be non-negative.
    uint16_t num_of_samples = 100;

    try{
        // 1\ Calculation of the Linear Bézier curve p(t).
        std::vector<std::vector<float>> result_lin = Bezier::Linear<float, uint16_t>(points_lin, num_of_samples);
        // 2\ Calculation of the Quadratic Bézier Curve p(t).
        std::vector<std::vector<float>> result_quadratic = Bezier::Quadratic<float, uint16_t>(points_quadratic, num_of_samples);
        // 3\ Calculation of the the Cubic Bézier Curve p(t).
        std::vector<std::vector<float>> result_cubic = Bezier::Cubic<float, uint16_t>(points_cubic, num_of_samples);

        /**
         * 4\ Calculation of the N-Degree Bézier curve p(t).
         * 
         * Note: The result of the calculation depends on the simplification factor.
         */
        Bezier::N_Degree<float, uint16_t> Bezier_Ndeg{num_of_samples};
        std::vector<std::vector<float>> result_n_deg = Bezier_Ndeg.Solve(points_n_deg, 1);

        if(display_curve){
            std::vector<std::vector<float>> result = result_n_deg;

            // Display the result of the calculation
            std::cout << "[";
            for(const auto & row : result | boost::adaptors::indexed(0)){
                std::cout << "[";
                for(const auto & col : row.value() | boost::adaptors::indexed(0)){
                    std::cout << col.value();
                    if(col.index() != result[0].size() - 1){
                        std::cout << ",";
                    }
                }
                if(row.index() != result.size() - 1){
                    std::cout << "],";
                }else{
                    std::cout << "]";
                }
            }
            std::cout << "]" << std::endl;
        }
    }catch (std::invalid_argument& e){
        std::cerr << e.what() << std::endl;
        return -1;
    }

    return 0;
}