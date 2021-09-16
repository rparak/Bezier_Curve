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
 * File Name: Utils.h
 * Github   : https://github.com/rparak
 */

#include <iostream>
#include <vector>
#include <boost/range/adaptor/indexed.hpp>

namespace Utils{
    /*
    Description:
        Simple utility functions.
    */

    template <typename T>
    std::vector<std::vector<T>> Transpose(const std::vector<std::vector<T>> data) {
        /*  
        Description:
            The function changes the row elements into column elements and the column elements into row elements.
        
        Args:
            (1) data [std::vector<std::vector<T>>]: Input array.

        Returns:
            (1) parameter [[std::vector<std::vector<T>>]: Transpose (rotate) data from rows to columns.
        */

        // Allocates space for the result
        std::vector<std::vector<T>> result(data[0].size(), std::vector<T>(data.size()));

        for(const auto & row : data | boost::adaptors::indexed(0)){
            for(const auto & col : row.value()  | boost::adaptors::indexed(0)){
                result[col.index()][row.index()] = data[row.index()][col.index()];
            }
        }

        return result;
    }

    template <typename T, typename U>
    std::vector<T> Linspace(const T t_start, const T t_stop, const U num_of_samples) {
        /*  
        Description:
            The function returns evenly spaced numbers over a specified interval.
        
        Args:
            (1) t_start [T]: The starting value of the sequenc.
            (2) t_stop [T]: The end value of the sequence.
            (3) number_of_samples [int]: Number of samples to generate. Must be non-negative.

        Returns:
            (1) parameter [std::vector<T>]: Returns data evenly spaced samples, calculated over the interval [t_start, t_stop].
        */

        if(num_of_samples < 0) { throw std::invalid_argument("The number of samples must not be negative."); }

        // Allocates space for the result
        std::vector<T> result(num_of_samples);

        if(num_of_samples < 2){
            result[0] = t_stop;
            return result;
        }

        // Delta Time
        T dt = (t_stop - t_start) / static_cast<T>(num_of_samples - 1);

        typename std::vector<T>::iterator iter;
        T value;
        for (iter = result.begin(), value = t_start; iter != result.end(); ++iter, value += dt){
            *iter = value;
        }

        return result;
    }

    uint16_t Binomial_Coefficient(const uint16_t n, const uint16_t k){
        /*  
        Description:
                Calculation binomial coofecient C, from pair of integers n ≥ k ≥ 0 and is written (n k). The binomial coefficients are the positive integers that occur as coefficients in the binomial theorem.
                
                (n k) = n! / (k! * (n - k)!)
                ...
                Simplification of the calculation:
                (n k) = ((n - k + 1) * (n - k + 2) * ... * (n - 1) * (n)) / (1 * 2 * ... * (k - 1) * k)
        
        Args:
            (1) n [uint16_t]: Integer number 1 (numerator)
            (2) k [uint16_t]: Integer number 2 (denumerator)

        Returns:
            (1) parameter [uint16_t]: Binomial coofecient C(n k).
        */

        if(k > n) { throw std::invalid_argument("The number n must be larger than or equal to k."); }

        if(k == 0){
            return 1;
        }else if(k == 1){
            return n;
        }else{
            uint16_t c_nk = 1;
            
            for(auto i = 0; i < k; ++i){
                c_nk *= (n - i);
                c_nk /= (i + 1);
            }

            return c_nk;
        }
    }
}
