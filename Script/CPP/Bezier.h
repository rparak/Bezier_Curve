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
 * File Name: Bezier.h
 * Github   : https://github.com/rparak
 */

#include "Utils.h"
#include <cmath>

// Initialization of constants
#define NUM_OF_ENTRY_POINTS_LINEAR    2
#define NUM_OF_ENTRY_POINTS_QUADRATIC 3
#define NUM_OF_ENTRY_POINTS_CUBIC     4
// Time t ∈ [0: The starting value of the sequence, 
//           1: The end value of the sequence]
#define T_START 0
#define T_STOP  1

namespace Bezier{
    /*
    Description:
        A Bézier curve is a parametric curve used in computer graphics and related fields.
    */

    template <typename T, typename U>
    std::vector<std::vector<T>> Linear(const std::vector<std::vector<T>> points, const U num_of_samples){
        /*
        Description:
            Given two control points p_{0} and p_{1} we define the linear Bézier curve to be the curve parametrized by:

            p(t) = (1 - t)*p_{0} + t*p_{1}, t ∈ [0, 1]
            
        Args:
            (1) num_of_samples [U]: Number of samples to generate. Must be non-negative.
            (2) points [p_{0, 1}] [std::vector<std::vector<T>>]: Multiple points to create a curve.

        Returns:
            (1) parameter [{0 .. Number of dimensions - 1}] [std::vector<std::vector<T>>]: Resulting points of the curve.

        Example:
            std::vector<std::vector<float>> result = Bezier::Linear<float>(points, number_of_samples)),
                where points are equal to [[px_id_0, py_id_0], [px_id_1, py_id_1]] in 2D space 
                and [[px_id_0, py_id_0, pz_id_0], [px_id_1, py_id_1, pz_id_1]] in 3D space
	    */
        if(points.size() != NUM_OF_ENTRY_POINTS_LINEAR) { throw std::invalid_argument("Insufficient number of entry points."); }

        // Return evenly spaced numbers over a specified interval
        std::vector<T> time = Utils::Linspace<T, U>(T_START, T_STOP, num_of_samples);

        // Transpose (rotate) data from rows to columns
        std::vector<std::vector<T>> points_T = Utils::Transpose(points);

        // Allocates space for the result
        std::vector<std::vector<T>> result(time.size(), std::vector<T>(points_T.size()));

        for(const auto & p : points_T | boost::adaptors::indexed(0)){
            for(const auto & t : time | boost::adaptors::indexed(0)){
                result[t.index()][p.index()] = (1 - t.value()) * p.value()[0] + t.value() * p.value()[1];
            }
        }

        return result;
    }

    template <typename T, typename U>
    std::vector<std::vector<T>> Quadratic(const std::vector<std::vector<T>> points, const U num_of_samples){
        /*
        Description:
            Given three control points p_{0}, p_{1} and p_{2} we define the quadratic Bézier curve (degree 2 Bézier curve)
            to be the curve parametrized by:

            p(t) = ((1 - t)^2)*p_{0} + 2*t*(1 - t)*p_{1} + (t^2)*p_{2}, t ∈ [0, 1]
            
        Args:
            (1) num_of_samples [U]: Number of samples to generate. Must be non-negative.
            (2) points [p_{0, 1, 3}] [std::vector<std::vector<T>>]: Multiple points to create a curve.

        Returns:
            (1) parameter [{0 .. Number of dimensions - 1}] [std::vector<std::vector<T>>]: Resulting points of the curve.

        Example:
            std::vector<std::vector<float>> result = Bezier::Quadratic<float>(points, number_of_samples)),
                where points are equal to [[px_id_0, py_id_0], [px_id_1, py_id_1], [px_id_2, py_id_2]] in 2D space 
                and [[px_id_0, py_id_0, pz_id_0], [px_id_1, py_id_1, pz_id_1], [px_id_2, py_id_2, pz_id_2]] in 3D space
	    */

        if(points.size() != NUM_OF_ENTRY_POINTS_QUADRATIC) { throw std::invalid_argument("Insufficient number of entry points."); }

        // Return evenly spaced numbers over a specified interval
        std::vector<T> time = Utils::Linspace<T, U>(T_START, T_STOP, num_of_samples);

        // Transpose (rotate) data from rows to columns
        std::vector<std::vector<T>> points_T = Utils::Transpose(points);

        // Allocates space for the result
        std::vector<std::vector<T>> result(time.size(), std::vector<T>(points_T.size()));

        for(const auto & p : points_T | boost::adaptors::indexed(0)){
            for(const auto & t : time | boost::adaptors::indexed(0)){
                result[t.index()][p.index()] = ((1 - t.value()) * (1 - t.value())) * p.value()[0] + 2 * t.value() * (1 - t.value()) * p.value()[1] +  (t.value() * t.value()) * p.value()[2];
            }
        }

        return result;
    }

    template <typename T, typename U>
    std::vector<std::vector<T>> Cubic(const std::vector<std::vector<T>> points, const U num_of_samples){
        /*
        Description:
            Given four control points p_{0}, p_{1}, p_{2} and p_{3} we define the cubic Bézier curve (degree 3 Bézier curve) to
            be the curve parametrized by:

            p(t) = ((1 - t)^3)*p_{0} + 3*t*((1 - t)^2)*p_{1} + (3*t^2)*(1 - t)*p_{2} + (t^3) * p_{3}, t ∈ [0, 1]
            
        Args:
            (1) num_of_samples [U]: Number of samples to generate. Must be non-negative.
            (2) points [p_{0, 1, 2, 3}] [std::vector<std::vector<T>>]: Multiple points to create a curve.

        Returns:
            (1) parameter [{0 .. Number of dimensions - 1}] [std::vector<std::vector<T>>]: Resulting points of the curve.

        Example:
            std::vector<std::vector<float>> result = Bezier::Cubic<float>(points, number_of_samples)),
                where points are equal to [[px_id_0, py_id_0], [px_id_1, py_id_1], [px_id_2, py_id_2], [px_id_3, py_id_3]] in 2D space 
                and [[px_id_0, py_id_0, pz_id_0], [px_id_1, py_id_1, pz_id_1], [px_id_2, py_id_2, pz_id_2], [px_id_3, py_id_3, pz_id_2]] in 3D space
	    */

        if(points.size() != NUM_OF_ENTRY_POINTS_CUBIC) { throw std::invalid_argument("Insufficient number of entry points."); }

        // Return evenly spaced numbers over a specified interval
        std::vector<T> time = Utils::Linspace<T, U>(T_START, T_STOP, num_of_samples);

        // Transpose (rotate) data from rows to columns
        std::vector<std::vector<T>> points_T = Utils::Transpose(points);

        // Allocates space for the result
        std::vector<std::vector<T>> result(time.size(), std::vector<T>(points_T.size()));

        for(const auto & p : points_T | boost::adaptors::indexed(0)){
            for(const auto & t : time | boost::adaptors::indexed(0)){
                result[t.index()][p.index()] = ((1 - t.value()) * (1 - t.value()) * (1 - t.value())) * p.value()[0] + 3 * t.value() * ((1 - t.value()) * (1 - t.value())) * p.value()[1] + 3 * (t.value() * t.value()) * (1 - t.value()) * p.value()[2] + (t.value() * t.value() * t.value()) * p.value()[3];
            }
        }

        return result;
    }

    template <typename T, typename U>
    class N_Degree {
        /*
        Description:
            Class for efficient solution of N-degree Bézier curve.
            
        Initialization of the Class:
            Input:
                (1) num_of_samples [U]: Number of samples to generate. Must be non-negative.

        Example:
            Initialization:
                Bezier::N_Degree<float, int> Cls{num_of_samples};
            Calculation:
                std::vector<std::vector<float>> result = Cls.Solve(points, simplification_factor);

                where points are equal to [[px_id_0, py_id_0], .., [px_id_n, py_id_n]] in 2D space 
                and [[px_id_0, py_id_0, pz_id_0], .., [px_id_n, py_id_n, pz_id_n]] in 3D space
        */

        public:
            // Initialization of input parameters
            explicit N_Degree(const U num_of_samples) 
                : __num_of_samples(num_of_samples) 
            {}

            // Function for automatic calculation of a suitably selected Bézier curve
            std::vector<std::vector<T>> Solve(const std::vector<std::vector<T>> points, uint16_t simplification_factor);
        private:
            // Number of samples to generate. Must be non-negative.
            U __num_of_samples{0};
            // Mumber of points in the matrix: 2D - 2 points, 3D - 3 points
            U __number_of_points;
            // Binomial Coefficient
            uint16_t __c_nk;

            // Return evenly spaced numbers over a specified interval.
            std::vector<T> __time = Utils::Linspace<T, U>(T_START, T_STOP, this->__num_of_samples);
            // Function to simplify the path through the simplification factor.
            std::vector<std::vector<T>> __Path_Simplification(const std::vector<std::vector<T>> points, uint16_t simplification_factor);
            // Function to calculate position in a specific iteration.
            std::vector<std::vector<T>> __N_Index(const std::vector<T> point, const U index, const uint16_t c_ni);
            // The main control function for creating a Bézier curve of degree n.
            std::vector<std::vector<T>> __N_Degree(const std::vector<std::vector<T>> points);
    };

    template <typename T, typename U>
    std::vector<std::vector<T>> N_Degree<T, U>::__Path_Simplification(const std::vector<std::vector<T>> points, uint16_t simplification_factor){
        /*
        Description:
            Function to simplify the path through the simplification factor. The first and end points do not change, the others 
            depend on the factor coefficient.

            Example:
                Input Points: 
                    points = [1.0, 1.0], [1.25, 2.0], [1.75, 2.0], [2.0, 1.0], [1.0, -1.0], [1.25, -2.0], [1.75, -2.0], [2.0, -1.0]
                Number of points: 
                    n = 8
                Simplification Factor:
                    1\ Example:
                        simplification_factor = 1
                        points_new = [1.0, 1.0], [1.25, 2.0], [1.75, 2.0], [2.0, 1.0], [1.0, -1.0], [1.25, -2.0], [1.75, -2.0], [2.0, -1.0]
                        n = 8

                    2\ Example:
                        simplification_factor = 2
                        points_new = [1.0, 1.0], [None], [1.75, 2.0], [None], [1.0, -1.0], [None], [1.75, -2.0], [2.0, -1.0] 
                        points_new = [1.0, 1.0], [1.75, 2.0], [1.0, -1.0], [1.75, -2.0], [2.0, -1.0]
                        n = 5

        Args:
            (1) points [p_{0, .., n}] [std::vector<std::vector<T>>]: Multiple points to create a curve.
            (2) simplification_factor [uint16_t]: Simplification factor for the simplify the path.

        Return:
            (1) parameter [{0 .. Number of dimensions - 1}] [std::vector<std::vector<T>>]: New simplified matrix of points to create a curve.

        */

        // Allocates space for the result
        std::vector<std::vector<T>> result{0};

        result.push_back(points[0]);

        for(auto i = 1; i < points.size() - 1; ++i){
            if(i % simplification_factor == 0){
               result.push_back(points[i]); 
            }
        }

        if(result.back() != points.back()){ 
            result.push_back(points.back());
        }

        return result;
    }

    template <typename T, typename U>
    std::vector<std::vector<T>> N_Degree<T, U>::__N_Index(const std::vector<T> point, const U i, const uint16_t c_ni){
        /*
        Description: 
            Given n + 1 control points p_{0}, p_{1},..., p_{n} we define the degree n Bezier curve to
            be the curve parametrized by (De Casteljau's algorithm):

            p(t) = sum(i = 0 -> n) (C(n i)) * (t ^ i) * ((1 - t) ^ (n - i)) * p_{i}, t ∈ [0, 1]

            where C(n i) is a binomial coefficient.

        Args:
            (1) point [std::vector<T>]: Point (2D/3D) in interation (i).
            (2) i [U]: Iteration.
            (3) c_ni [uint16_t]: Binomial coofecient C(n i) in iteration (i).

        Returns:
            (1) parameter [{0 .. Number of dimensions - 1}] [std::vector<std::vector<T>>]: Resulting points of the curve.
        */

        // Allocates space for the result
        std::vector<std::vector<T>> result(this->__time.size(), std::vector<T>(point.size()));

        for(const auto & p : point | boost::adaptors::indexed(0)){
            for(const auto & t : this->__time | boost::adaptors::indexed(0)){
                result[t.index()][p.index()] = c_ni * std::pow(t.value(), i) * std::pow(1 - t.value(), this->__number_of_points - i) * p.value();
            }
        }

        return result;
    }

    template <typename T, typename U>
    std::vector<std::vector<T>> N_Degree<T, U>::__N_Degree(const std::vector<std::vector<T>> points){
        /*
        Description: 
            The main control function for creating a Bézier curve of degree n.

        Args:
            (1) points [p_{0, .., n}] [std::vector<std::vector<T>>]: Multiple points to create a curve.

        Returns:
            (1) parameter [{0 .. Number of dimensions - 1}] [std::vector<std::vector<T>>]: Resulting points of the curve.
        */

        // Allocates space for the variables
        std::vector<std::vector<T>> result(this->__time.size(), std::vector<T>(points.size()));
        std::vector<std::vector<T>> aux_result{0};

        // Mumber of points in the matrix
        this->__number_of_points = points.size() - 1;

        // Calculation of the first iteration
        uint16_t c_nk = Utils::Binomial_Coefficient(this->__number_of_points, 0);
        result = this->__N_Index(points[0], 0, c_nk);

        for(auto i = 1; i < this->__number_of_points + 1; ++i){
            // Binomial cooficient in interation (i)
            c_nk = Utils::Binomial_Coefficient(this->__number_of_points, i);

            // Calculation positions in iteration (i)
            aux_result = this->__N_Index(points[i], i, c_nk);

            // The sum of all positions for the resulting Bézier curve
            for(const auto & row : aux_result | boost::adaptors::indexed(0)){
                for(const auto & col : row.value()  | boost::adaptors::indexed(0)){
                    result[row.index()][col.index()] += aux_result[row.index()][col.index()];
                }
            }
        }

        return result;
    }

    template <typename T, typename U>
    std::vector<std::vector<T>> N_Degree<T, U>::Solve(const std::vector<std::vector<T>> points, uint16_t simplification_factor) {
        /*
        Description:
            Function for automatic calculation of a suitably selected Bézier curve.
        
        Args:
            (1) points [p_{0, .., n}] [std::vector<std::vector<T>>]: Multiple points to create a curve.
            (2) simplification_factor [uint16_t]: Simplification factor for the simplify the path.

        Returns:
            (1) parameter [{0 .. Number of dimensions - 1}] [std::vector<std::vector<T>>]: Resulting points of the curve.
        */

        if(points.size() < NUM_OF_ENTRY_POINTS_LINEAR) { throw std::invalid_argument("Insufficient number of entry points."); }

        // Allocates space for the points (new)
        std::vector<std::vector<T>> points_new(this->__time.size(), std::vector<T>(points.size()));

        /*
        Description:
            If the number of input points is greater than 4 and variable simplification_factor > 1, the program chooses the n_points calculation method. But if the simplification 
            coefficient is greater than or equal to 1, the program can choose another method and the principle of calculation will be faster.
        */
        if(simplification_factor > 1 && points.size() > 4){
            // If the coefficient coefficient is greater than 1, simplify the path
            points_new = this->__Path_Simplification(points, simplification_factor);
        }else{
            points_new = points;
        }

        /*
        Description:
            Selects the calculation method based on the number of points in the matrix (p).
        */
        switch(points_new.size()){
            case 2:
            {
                return Bezier::Linear<T, U>(points_new, this->__num_of_samples);
            }
            break;

            case 3:
            {
                return Bezier::Quadratic<T, U>(points_new, this->__num_of_samples);
            }
            break;

            case 4:
            {
                return Bezier::Cubic<T, U>(points_new, this->__num_of_samples);
            }
            break;

            default:
            {
                return this->__N_Degree(points_new);
            }
            break;
        }
    }
}