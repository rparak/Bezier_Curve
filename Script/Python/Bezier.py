"""
## =========================================================================== ## 
MIT License
Copyright (c) 2021 Roman Parak
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
## =========================================================================== ## 
Author   : Roman Parak
Email    : Roman.Parak@outlook.com
Github   : https://github.com/rparak
File Name: Bezier.py
## =========================================================================== ## 
"""

# Numpy (Array computing Lib.) [pip3 install numpy]
import numpy as np

# Support for type hints
import typing

# Initialization of constants
CONST_NUM_OF_ENTRY_POINTS_LINEAR    = 2
CONST_NUM_OF_ENTRY_POINTS_QUADRATIC = 3
CONST_NUM_OF_ENTRY_POINTS_CUBIC     = 4
# Time t ∈ [0: The starting value of the sequence, 
#           1: The end value of the sequence]
CONST_T_START = 0
CONST_T_STOP  = 1

def Linear(num_of_samples: typing.Union[int], points: typing.Union[typing.List[int], typing.List[float]]) -> typing.Union[typing.List[int], typing.List[float]]: 
    """
    Description:
        Given two control points p_{0} and p_{1} we define the linear Bézier curve to be the curve parametrized by:

        p(t) = (1 - t)*p_{0} + t*p_{1}, t ∈ [0, 1]

    Args:
        (1) num_of_samples [INT]: Number of samples to generate. Must be non-negative.
        (2) points [p_{0, 1}] [Int/Float Matrix]: Multiple points to create a curve.

    Returns:
        (1) parameter [{0 .. Number of dimensions - 1}] [Int/Float Matrix]: Resulting points of the curve.

    Example:
        res = Linear(num_of_samples, points),
            where points are equal to [[px_id_0, py_id_0], [px_id_1, py_id_1]] in 2D space 
            and [[px_id_0, py_id_0, pz_id_0], [px_id_1, py_id_1, pz_id_1]] in 3D space
    """
    try:
        assert len(points) == CONST_NUM_OF_ENTRY_POINTS_LINEAR
        assert(num_of_samples >= 0)

        # Return evenly spaced numbers over a specified interval.
        t = np.linspace(CONST_T_START, CONST_T_STOP, num_of_samples)

        return [(1 - t) * p[0] + t * p[1] 
                for _, p in enumerate(np.transpose(points))]

    except AssertionError as error:
        print('[ERROR] Insufficient number of entry points.')
        print('[ERROR] The correct number of entry points is %d.' % CONST_NUM_OF_ENTRY_POINTS_LINEAR)
        print('[ERROR] The number of samples must not be negative.')

def Quadratic(num_of_samples: typing.Union[int], points: typing.Union[typing.List[int], typing.List[float]]) -> typing.Union[typing.List[int], typing.List[float]]: 
    """
    Description:
        Given three control points p_{0}, p_{1} and p_{2} we define the quadratic Bézier curve (degree 2 Bézier curve)
        to be the curve parametrized by:

        p(t) = ((1 - t)^2)*p_{0} + 2*t*(1 - t)*p_{1} + (t^2)*p_{2}, t ∈ [0, 1]

    Args:
        (1) num_of_samples [INT]: Number of samples to generate. Must be non-negative.
        (2) points [p_{0, 1, 2}] [Int/Float Matrix]: Multiple points to create a curve.

    Returns:
         (1) parameter [{0 .. Number of dimensions - 1}] [Int/Float Matrix]: Resulting points of the curve.

    Example:
        res = Quadratic(t, p),
            where points are equal to [[px_id_0, py_id_0], [px_id_1, py_id_1], [px_id_2, py_id_2]] in 2D space 
            and [[px_id_0, py_id_0, pz_id_0], [px_id_1, py_id_1, pz_id_1], [px_id_2, py_id_2, pz_id_2]] in 3D space
    """

    try:
        assert len(points) == CONST_NUM_OF_ENTRY_POINTS_QUADRATIC
        assert(num_of_samples >= 0)

        # Return evenly spaced numbers over a specified interval.
        t = np.linspace(CONST_T_START, CONST_T_STOP, num_of_samples)

        return [(1 - t)**2 * p[0] + 2 * t * (1 - t) * p[1] + t**2 * p[2] 
                for _, p in enumerate(np.transpose(points))]

    except AssertionError as error:
        print('[ERROR] Insufficient number of entry points.')
        print('[ERROR] The correct number of entry points is %d.' % CONST_NUM_OF_ENTRY_POINTS_QUADRATIC)
        print('[ERROR] The number of samples must not be negative.')

def Cubic(num_of_samples: typing.Union[int], points: typing.Union[typing.List[int], typing.List[float]]) -> typing.Union[typing.List[int], typing.List[float]]: 
    """
    Description:
        Given four control points p_{0}, p_{1}, p_{2} and p_{3} we define the cubic Bézier curve (degree 3 Bézier curve) to
        be the curve parametrized by:

        p(t) = ((1 - t)^3)*p_{0} + 3*t*((1 - t)^2)*p_{1} + (3*t^2)*(1 - t)*p_{2} + (t^3) * p_{3}, t ∈ [0, 1]

    Args:
        (1) num_of_samples [INT]: Number of samples to generate. Must be non-negative.
        (2) points [p_{0, 1, 2, 3}] [Int/Float Matrix]: Multiple points to create a curve.

    Returns:
        (1) parameter [{0 .. Number of dimensions - 1}] [Int/Float Matrix]: Resulting points of the curve.

    Example:
        res = Cubic(t, p),
            where points are equal to [[px_id_0, py_id_0], [px_id_1, py_id_1], [px_id_2, py_id_2], [px_id_3, py_id_3]] in 2D space 
            and [[px_id_0, py_id_0, pz_id_0], [px_id_1, py_id_1, pz_id_1], [px_id_2, py_id_2, pz_id_2], [px_id_3, py_id_3, pz_id_2]] in 3D space
    """
    try:
        assert len(points) == CONST_NUM_OF_ENTRY_POINTS_CUBIC
        assert(num_of_samples >= 0)

        # Return evenly spaced numbers over a specified interval.
        t = np.linspace(CONST_T_START, CONST_T_STOP, num_of_samples)

        return [((1 - t)**3) * (p[0]) + (3 * t * (1 - t)**2) * (p[1]) + 3 * (t**2) * (1 - t) * p[2] + (t**3) * p[3] 
                for _, p in enumerate(np.transpose(points))]

    except AssertionError as error:
        print('[ERROR] Insufficient number of entry points.')
        print('[ERROR] The correct number of entry points is %d.' % CONST_NUM_OF_ENTRY_POINTS_CUBIC)
        print('[ERROR] The number of samples must not be negative.')

class N_Degree(object):
    """
    Description:
        Class for efficient solution of N-degree Bézier curve.

        Note:
            A Bézier curve is a parametric curve used in computer graphics and related fields.

    Initialization of the Class:
        Input:
            (1) num_of_samples [INT]: Number of samples to generate. Must be non-negative.

        Example:
            Initialization:
                Cls = N_Degree(num_of_samples)
            Calculation:
                res = Cls.Solve(p, simplification_factor)
            
                where p is equal to [[px_id_0, py_id_0], .., [px_id_n, py_id_n]] in 2D space 
                and [[px_id_0, py_id_0, pz_id_0], .., [px_id_n, py_id_n, pz_id_n]] in 3D space
    """

    def __init__(self, num_of_samples: typing.Union[int]) -> None:
        # << PUBLIC >> #
        try:
            assert(num_of_samples >= 0)

            # Return evenly spaced numbers over a specified interval.
            self.t = np.linspace(CONST_T_START, CONST_T_STOP, num_of_samples)

        except AssertionError as error:
            print('[ERROR] The number of samples must not be negative.')

        # << PRIVATE >> #
        # Points [Float Matrix]
        self.__points = []
        # Number of samples to generate
        self.__num_of_samples = num_of_samples

    @staticmethod
    def __path_simplification(points, simplification_factor):
        """
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
            (1) simplification_factor [INT]: Simplification factor for the simplify the path.

        Return:
            (1) parameter{1} [Int/Float Matrix]: New simplified matrix of points to create a curve.

        """

        points_aux = []

        points_aux.append(points[0])

        for i in range(1, len(p) - 1):
            if i % simplification_factor == 0:
                points_aux.append(p[i])

        if points_aux[len(points_aux) - 1] != p[len(p) - 1]:
            points_aux.append(points[len(points) - 1])

        return points_aux
    
    @staticmethod
    def __binomial_coefficient(n, k):
        """
        Description:
            Calculation binomial coofecient C, from pair of integers n ≥ k ≥ 0 and is written (n k). The binomial coefficients are the positive integers that occur as coefficients in the binomial theorem.

            (n k) = n! / (k! * (n - k)!)
            ...
            Simplification of the calculation:
            (n k) = ((n - k + 1) * (n - k + 2) * ... * (n - 1) * (n)) / (1 * 2 * ... * (k - 1) * k)
        
        Args:
            (1) n [INT]: Integer number 1 (numerator)
            (2) k [INT]: Integer number 2 (denumerator)

        Returns:
            (1) parameter{1} [INT]: Binomial coofecient C(n k).
        """
        try:
            assert(n >= k)
            
            c_nk = 1
            
            if k == 0 or k == 1:
                return c_nk;
            else:
                for i in range(0, k):
                    c_nk *= (n - i)
                    c_nk /= (i + 1)
                    
                return c_nk

        except AssertionError as error:
            print('[ERROR] The number n must be larger than or equal to k.')
            return 0

    def __n_index_curve(self, idx, point, n, c_ni):
        """
        Description:
            Given n + 1 control points p_{0}, p_{1},..., p_{n} we define the degree n Bezier curve to
            be the curve parametrized by (De Casteljau's algorithm):

            p(t) = sum(i = 0 -> n) (C(n i)) * (t ^ i) * ((1 - t) ^ (n - i)) * p_{i}, t ∈ [0, 1]

            where C(n i) is a binomial coefficient.

        Args:
            (1) idx [INT]: Iteration.
            (2) point [Int/Float Matrix]: Point (2D/3D) in interation (i).
            (3) n [INT]: Number of points.
            (4) c_ni [INT]: Binomial coofecient C(n i) in iteration (i).

        Returns:
            (1) parameter{1 .. self.__num_of_dimensions} [Int/Float Matrix]: Results of curve values in iteration (i).
        """

        return [c_ni * (self.t**idx) * ((1 - self.t)**(n - idx)) * p 
                for _, p in enumerate(point)]

    def __n_degree(self):
        """
        Description: 
            The main control function for creating a Bézier curve of degree n.

        Returns:
            (1) parameter [{0 .. Number of dimensions - 1}] [Int/Float Matrix]: Resulting points of the curve.
        """
        
        # Number of points in matrix
        n = len(self.__points) - 1

        # Calculation of binomial cooficient of the first iteration
        c_nk = self.__binomial_coefficient(n, 0)

        # Calculation of the first curve positions
        result = self.__n_index_curve(0, self.__points[0], n, c_nk)

        for i in range(1, n + 1):
            # Binomial cooficient in interation (i)
            c_nk = self.__binomial_coefficient(n, i)

            # Calculation positions in iteration (i)
            aux_result = self.__n_index_curve(i, self.__points[i], n, c_nk)

            # The sum of all positions for the resulting Bézier curve
            for j in range(0, len(aux_result)):
                result[j] += aux_result[j]
            
        return result

    def Solve(self, points: typing.Union[typing.List[int], typing.List[float]], simplification_factor: typing.Union[int]) -> typing.Union[typing.List[int], typing.List[float]]:
        """
        Description:
            Function for automatic calculation of a suitably selected Bézier curve.

        Args:
            (1) points [p_{0, .., n}] [Int/Float Matrix]: Multiple points to create a curve.
            (2) simplification_factor [INT]: Simplification factor for the simplify the path.

        Return:
            (1) parameter [{0 .. Number of dimensions - 1}] [Int/Float Matrix]: Resulting points of the curve.
        """

        try:
            assert len(points) > 1

            # If the number of input points is greater than 4 and variable simplification_factor > 1, the program chooses the n_points calculation method. But if the simplification 
            # coefficient is greater than or equal to 1, the program can choose another method and the principle of calculation will be faster.
            if simplification_factor > 1 and len(points) > 4:
                # If the coefficient coefficient is greater than 1, simplify the path
                self.__points = self.__path_simplification(points, simplification_factor)
            else:
                self.__points = points

            # Selects the calculation method based on the number of points in the matrix (p).
            if len(self.__points) > 4:
                return self.__n_degree()
            if len(self.__points) == 4:
                return Cubic(self.__num_of_samples, self.__points)
            elif len(self.__points) == 3:
                return Quadratic(self.__num_of_samples, self.__points)
            elif len(self.__points) == 2:
                return Linear(self.__num_of_samples, self.__points)

        except AssertionError as error:
            print('[ERROR] Insufficient number of entry points.')
            print('[ERROR] The minimum number of entry points is %d.' % CONST_NUM_OF_ENTRY_POINTS_LINEAR)
