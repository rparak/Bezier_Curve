"""
## =========================================================================== ## 
MIT License
Copyright (c) 2020 Roman Parak
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
File Name: bezier_curve.py
## =========================================================================== ## 
"""

# Numpy (Array computing Lib.) [pip3 install numpy]
import numpy as np

class bezier_ctrl(object):
    """
    Description:
        A Bézier curve is a parametric curve used in computer graphics and related fields.

        The class shows several types of Bézier curves (Linear, Quadratic, Cubic, N-Degree).

    Initialization of the Class:

        Input:
            (1) Time Step  [INT]

        Example:
            Initialization:
                cls = bezier_ctrl(100)
            Calculation:
                res = cls.linear_curve(points)
                res = cls.quadratic_curve(points)
                res = cls.cubic_curve(points)
                res = cls.auto(points, simplification_factor)
            
            where the points are equal to [[1.0, 1.0], [1.25, 2.0]] in 2D space and [[1.0, 1.0], [1.25, 2.0], [1.75, 2.0]] in 3D space
    """

    def __init__(self, step):
        # << PUBLIC >> #
        # Time t ∈ [0, 1]
        self.t = np.linspace(0.0, 1.0, step)
        # << PRIVATE >> #
        # Points [Float Array]
        self.__p = []

    def linear_curve(self, p):
        """
        Description:
            Given two control points p_{0} and p_{1} we define the linear Bezier curve to be the curve parametrized by:

            p(t) = (1 - t)*p_{0} + t*p_{1}, t ∈ [0, 1]

        Args:
            (1) points [p_{0, 1}] [Float Array]: Multiple points to create a curve.

        Returns:
            (1) parameter{0 .. Number of dimensions - 1} [Float Array]: Results of curve values.
        """

        return [(1 - self.t) * point[0] + self.t * point[1] 
                for _, point in enumerate(np.transpose(p))]

    def quadratic_curve(self, p):
        """
        Description:
            Given three control points p_{0}, p_{1} and p_{2} we define the quadratic Bezier curve (degree 2 Bezier curve)
            to be the curve parametrized by:

            p(t) = ((1 - t)^2)*p_{0} + 2*t*(1 - t)*p_{1} + (t^2)*p_{2}, t ∈ [0, 1]

        Args:
            (1) points [p_{0, 1, 2}] [Float Array]: Multiple points to create a curve.

        Returns:
            (1) parameter{0 .. Number of dimensions - 1} [Float Array]: Results of curve values.
        """

        return [(1 - self.t)**2 * point[0] + 2 * self.t * (1 - self.t) * point[1] + self.t**2 * point[2] 
                for _, point in enumerate(np.transpose(p))]   

    def cubic_curve(self, p):
        """
        Description:
            Given four control points p_{0}, p_{1}, p_{2} and p_{3} we define the cubic Bezier curve (degree 3 Bezier curve) to
            be the curve parametrized by:

            p(t) = ((1 - t)^3)*p_{0} + 3*t*((1 - t)^2)*p_{1} + (3*t^2)*(1 - t)*p_{2} + (t^3) * p_{3}, t ∈ [0, 1]

        Args:
            (1) points [p_{0, 1, 2, 3}] [Float Array]: Multiple points to create a curve.

        Returns:
            (1) parameter{0 .. Number of dimensions - 1} [Float Array]: Results of curve values.
        """

        return [((1 - self.t)**3) * (point[0]) + (3 * self.t * (1 - self.t)**2) * (point[1]) + 3 * (self.t**2) * (1 - self.t) * point[2] + (self.t**3) * point[3] 
                for _, point in enumerate(np.transpose(p))]

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

        # Take advantage of symmetry
        if k > (n - k):
            k = n - k

        c_nk = 1

        # Calculation from the simplification equation
        for i in range(0, k):
            c_nk *= (n - i) # numerator
            c_nk /= (i + 1) # denumerator

        return c_nk

    def __n_index_curve(self, idx, p, n, c_ni):
        """
        Description:
            Given n + 1 control points p_{0}, p_{1},..., p_{n} we define the degree n Bezier curve to
            be the curve parametrized by (De Casteljau's algorithm):

            p(t) = sum(i = 0 -> n) (C(n i)) * (t ^ i) * ((1 - t) ^ (n - i)) * p_{i}, t ∈ [0, 1]

            where C(n i) is a binomial coefficient.

        Args:
            (1) idx [INT]: Iteration.
            (2) p [Float Array]: Point (2D/3D) in interation (i).
            (3) n [INT]: Number of points.
            (4) c_ni [INT]: Binomial coofecient C(n i) in iteration (i).

        Returns:
            (1) parameter{1 .. self.__num_of_dimensions} [Float Array]: Results of curve values in iteration (i).
        """

        return [c_ni * (self.t**idx) * ((1 - self.t)**(n - idx)) * point 
                for _, point in enumerate(p)]

    @staticmethod
    def __path_simplification(p, simplification_factor):
        """
        Description:
            Function to simplify the path through the simplification factor. The first and end points do not change, the others 
            depend on the factor coefficient.

            Example:
                Input Points: 
                    p = [1.0, 1.0], [1.25, 2.0], [1.75, 2.0], [2.0, 1.0], [1.0, -1.0], [1.25, -2.0], [1.75, -2.0], [2.0, -1.0]
                Number of points: 
                    n = 8
                Simplification Factor:
                    s_f           = 1
                    p_aux (new p) = [1.0, 1.0], [1.25, 2.0], [1.75, 2.0], [2.0, 1.0], [1.0, -1.0], [1.25, -2.0], [1.75, -2.0], [2.0, -1.0]
                    n             = 8

                    s_f           = 2
                    p_aux (new p) = [1.0, 1.0], [None], [1.75, 2.0], [None], [1.0, -1.0], [None], [1.75, -2.0], [2.0, -1.0] 
                    p_aux (new p) = [1.0, 1.0], [1.75, 2.0], [1.0, -1.0], [1.75, -2.0], [2.0, -1.0]
                    n             = 5
        Args:
            (1) simplification_factor [INT]: Simplification factor for the simplify the path.

        Return:
            (1) parameter{1} [Float Array]: New simplified array of points to create a curve.

        """

        p_aux = []

        p_aux.append(p[0])

        for i in range(1, len(p) - 1):
            if i % simplification_factor == 0:
                p_aux.append(p[i])

        if p_aux[len(p_aux) - 1] != p[len(p) - 1]:
            p_aux.append(p[len(p) - 1])

        return p_aux

    def __n_points(self):
        """
        Description: 
            The main control function for creating a Bézier curve of degree n.
        """
        
        # Number of points in array
        n = len(self.__p) - 1

        # Calculation of binomial cooficient of the first iteration
        c_nk = self.__binomial_coefficient(n, 0)

        # Calculation of the first curve positions
        result = self.__n_index_curve(0, self.__p[0], n, c_nk)

        for i in range(1, n + 1):
            # Binomial cooficient in interation (i)
            c_nk = self.__binomial_coefficient(n, i)

            # Calculation positions in iteration (i)
            aux_result = self.__n_index_curve(i, self.__p[i], n, c_nk)

            # The sum of all positions for the resulting Bézier curve
            for j in range(0, len(aux_result)):
                result[j] += aux_result[j]
            
        return result

    def auto(self, points, simplification_factor):
        """
        Description:
            Function of automatic calculation of individual Bézier curves. 
            
        Args:
            (1) points [p_{0, .., n}] [Float Array]: Multiple points to create a curve.
            (2) simplification_factor [INT]: Simplification factor for the simplify the path.

        Return:
            (1) parameter{0 .. Number of dimensions - 1} [Float Array]: Results of curve values.
        """

        try:
            assert len(points) > 1

            # If the number of input points is greater than 4 and variable simplification_factor > 1, the program chooses the n_points calculation method. But if the simplification 
            # coefficient is greater than or equal to 1, the program can choose another method and the principle of calculation will be faster.
            if simplification_factor > 1 and len(points) > 4:
                # If the coefficient coefficient is greater than 1, simplify the path
                self.__p = self.__path_simplification(points, simplification_factor)
            else:
                self.__p = points

            # Selects the calculation method based on the number of points in the array (p).
            if len(self.__p) > 4:
                return self.__n_points()
            if len(self.__p) == 4:
                return self.cubic_curve(self.__p)
            elif len(self.__p) == 3:
                return self.quadratic_curve(self.__p)
            elif len(self.__p) == 2:
                return self.linear_curve(self.__p)

        except AssertionError as error:
            print('[ERROR] Insufficient number of entry points.')
            print('[ERROR] The minimum number of entry points is 2.')
