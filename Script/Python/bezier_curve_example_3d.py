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
File Name: bezier_curve_example_3d.py
## =========================================================================== ## 
"""

# System (Default Lib.)
import sys
# Numpy (Array computing Lib.) [pip3 install numpy]
import numpy as np
# Mtaplotlib (Visualization Lib.) [pip3 install matplotlib]
import matplotlib.pyplot as plt


class bezier_ctrl(object):
    """
    Description:
        A Bézier curve is a parametric curve used in computer graphics and related fields.

        The class shows several types of Bézier curves (Linear, Quadratic, Cubic).
    """
    def __init__(self, p, step):
        # << PUBLIC >> #
        # Time t ∈ [0, 1]
        self.t = np.linspace(0.0, 1.0, step)
        # Points
        self.p = p
        # Path simplification factor (skip some points)
        self.simplification_factor = 1
        # << PRIVATE >> #
        # Display (Plot) variable.
        self.__plt  = plt
        # Axis (3D)
        self.__axis = self.__plt.axes(projection='3d')

    @staticmethod
    def __linear_curve(p_0, p_1, t):
        """
        Description:
            Given two control points p_{0} and p_{1} we define the linear Bezier curve to be the curve parametrized by:

            p(t) = (1 - t)*p_{0} + t*p_{1}, t ∈ [0, 1]

        Args:
            (1 - 2) p_0, p_0 [Float Array]: Multiple points to create a curve.
            (3) t [Float Array]: Time variable.

        Returns:
            (1 - 3) parameter{1}, parameter{2}, parameter{3} [Float Array]: Results of curve values.

        Examples:
            self.__linear_curve([1.0, 1.0], [2.0, 2.0])
        """

        x = (1 - t) * p_0[0] + t * p_1[0]
        y = (1 - t) * p_0[1] + t * p_1[1]
        z = (1 - t) * p_0[2] + t * p_1[2]

        return x, y, z

    @staticmethod
    def __quadratic_curve(p_0, p_1, p_2, t):
        """
        Description:
            Given three control points p_{0}, p_{1} and p_{2} we define the quadratic Bezier curve (degree 2 Bezier curve)
            to be the curve parametrized by:

            p(t) = ((1 - t)^2)*p_{0} + 2*t*(1 - t)*p_{1} + (t^2)*p_{2}, t ∈ [0, 1]

        Args:
            (1 - 2) p_0, p_0, p_2 [Float Array]: Multiple points to create a curve.
            (3) t [Float Array]: Time variable.

        Returns:
            (1 - 3) parameter{1}, parameter{2}, parameter{3} [Float Array]: Results of curve values.

        Examples:
            self.__quadratic_curve([1.0, 1.0], [2.0, 2.0], [3.0, 2.0])
        """

        x = (1 - t)**2 * p_0[0] + 2 * t * (1 - t) * p_1[0] + t**2 * p_2[0]
        y = (1 - t)**2 * p_0[1] + 2 * t * (1 - t) * p_1[1] + t**2 * p_2[1]
        z = (1 - t)**2 * p_0[2] + 2 * t * (1 - t) * p_1[2] + t**2 * p_2[2]

        return x, y, z   

    @staticmethod
    def __cubic_curve(p_0, p_1, p_2, p_3, t):
        """
        Description:
            Given four control points p_{0}, p_{1}, p_{2} and p_{3} we define the cubic Bezier curve (degree 3 Bezier curve) to
            be the curve parametrized by:

            p(t) = ((1 - t)^3)*p_{0} + 3*t*((1 - t)^2)*p_{1} + (3*t^2)*(1 - t)*p_{2} + (t^3) * p_{3}, t ∈ [0, 1]

        Args:
            (1 - 2) p_0, p_0, p_2 [Float Array]: Multiple points to create a curve.
            (3) t [Float Array]: Time variable.

        Returns:
            (1 - 3) parameter{1}, parameter{2}, parameter{3} [Float Array]: Results of curve values.

        Examples:
            self.__cubic_curve([1.0, 1.0], [2.0, 2.0], [3.0, 2.0], [4.0, 1.0])
        """

        x = ((1 - t)**3) * (p_0[0]) + (3 * t * (1 - t)**2) * (p_1[0]) + 3 * (t**2) * (1 - t) * p_2[0] + (t**3) * p_3[0]
        y = ((1 - t)**3) * (p_0[1]) + (3 * t * (1 - t)**2) * (p_1[1]) + 3 * (t**2) * (1 - t) * p_2[1] + (t**3) * p_3[1]
        z = ((1 - t)**3) * (p_0[2]) + (3 * t * (1 - t)**2) * (p_1[2]) + 3 * (t**2) * (1 - t) * p_2[2] + (t**3) * p_3[2]

        return x, y, z

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
        for i in range(k):
            c_nk *= (n - i) # numerator
            c_nk /= (i + 1) # denumerator

        return c_nk

    @staticmethod
    def __n_index_curve(i, p, t, n, c_ni):
        """
        Description:
            Given n + 1 control points p_{0}, p_{1},..., p_{n} we define the degree n Bezier curve to
            be the curve parametrized by (De Casteljau's algorithm):

            p(t) = sum(i = 0 -> n) (C(n i)) * (t ^ i) * ((1 - t) ^ (n - i)) * p_{i}, t ∈ [0, 1]

            where C(n i) is a binomial coefficient.

        Args:
            (1) i [INT]: Iteration.
            (2) p [Float Array]: Point (x, y) in interation (i).
            (3) t [Float Array]: Time variable.
            (4) n [INT]: Number of points.
            (5) c_ni [INT]: Binomial coofecient C(n i) in iteration (i).

        Returns:
            (1 - 2) parameter{1}, parameter{2}, parameter{3} [Float Array]: Results of curve values in iteration (i).
        """

        x = c_ni * (t**i) * ((1 - t)**(n - i)) * p[0]
        y = c_ni * (t**i) * ((1 - t)**(n - i)) * p[1]
        z = c_ni * (t**i) * ((1 - t)**(n - i)) * p[2]

        return x, y, z

    def __switch_dCtrl(self, s_index):
        """
        Description:
            Function to obtain a string with the number of points used in the calculation of the curve. (Figure Legend -> Label Name)

        Args:
            (1) s_index [INT]: Number of points for calculation.

        Returns:
            (1) parameter{1} [String]: The resulting string for the label.
        """

        # Switch Variable
        switch_var={
                2: r'Points: $p_{0}, p_{1}$',
                3: r'Points: $p_{0}, p_{1}, p_{2}$',
                4: r'Points: $p_{0}, p_{1}, p_{2}, p_{3}$',
        }

        # Return Result (Get the string with number of points for the Legend Label)
        return switch_var.get(s_index, r'Points: $p_{0}, p_{1},..., p_{n}$')

    def __two_points(self):
        """
        Description:
            Function to create a multiple linear Bézier curve from two points.
        """

        x, y, z = self.__linear_curve(self.p[0], self.p[1], self.t)

        # Display the Linear Bézier Curve p(t) -> x, y
        self.__axis.plot3D(x, y, z, 'r--', label=r'Linear Bezier Curve: [$p_{0}$, $p_{1}$]' , linewidth=2.5)

    def __three_points(self):
        """
        Description:
            Function to create multiple linear Bézier curves and a quadratic Bézier curve from three points.
        """

        for i in range(len(self.p) - 1):
            x, y, z = self.__linear_curve(self.p[i], self.p[i + 1], self.t)
            
            if i == (len(self.p) - 1) - 1:
                self.__axis.plot3D(x, y, z, 'r--', label=r'Linear Bezier Curve: [$p_{0}$, $p_{1}$]; [$p_{1}$, $p_{2}$]' , linewidth=2.5)
            else:
                self.__axis.plot3D(x, y, z, 'r--', linewidth=2.5)

        x, y, z = self.__quadratic_curve(self.p[0], self.p[1], self.p[2], self.t)

        self.__axis.plot3D(x, y, z, 'g--', label=r'Quadratic Bezier Curve: [$p_{0}$, $p_{1}$, $p_{2}$]', linewidth=2.5)

    def __four_points(self):
        """
        Description:
            Function to create multiple linear / quadratic Bézier curves and a cubic four-point Bézier curve.
        """

        for i in range(len(self.p) - 1):
            x, y, z  = self.__linear_curve(self.p[i], self.p[i + 1], self.t)
            
            if i == (len(self.p) - 1) - 1:
                self.__axis.plot3D(x, y, z, 'r--', label=r'Linear Bezier Curve: [$p_{0}$, $p_{1}$]; [$p_{1}$, $p_{2}$]; [$p_{2}$, $p_{3}$]' , linewidth=2.5)
            else:
                self.__axis.plot3D(x, y, z, 'r--', linewidth=2.5)

        for i in range(len(self.p) - 2):
            x, y, z = self.__quadratic_curve(self.p[i], self.p[i + 1], self.p[i + 2], self.t)

            if i == (len(self.p) - 2) - 1:
                self.__axis.plot3D(x, y, z, 'g--', label=r'Quadratic Bezier Curve: [$p_{0}$, $p_{1}$, $p_{2}$]; [$p_{1}$, $p_{2}$, $p_{3}$]', linewidth=2.5)
            else:
                self.__axis.plot3D(x, y, z, 'g--', linewidth=2.5)

        x, y, z = self.__cubic_curve(self.p[0], self.p[1], self.p[2], self.p[3], self.t)

        self.__axis.plot3D(x, y, z, 'b--', label=r'Cubic Bezier Curve: [$p_{0}$, $p_{1}$, $p_{2}$, $p_{3}$]', linewidth=2.5)

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
            Function to create multiple linear Bézier curves and a n-degree Bézier curve.
        """

        for i in range(len(self.p) - 1):
            x, y, z  = self.__linear_curve(self.p[i], self.p[i + 1], self.t)
            
            if i == (len(self.p) - 1) - 1:
                self.__axis.plot3D(x, y, z, 'r--', label=r'Linear Bezier Curve: [$p_{0}$, $p_{1}$]; [$p_{1}$, $p_{2}$]; ... ; [$p_{n - 1}$, $p_{n}$]' , linewidth=2.5)
            else:
                self.__axis.plot3D(x, y, z, 'r--', linewidth=2.5)

        # Number of points in array
        n = len(self.p) - 1

        # Calculation of binomial cooficient of the first iteration
        c_nk = self.__binomial_coefficient(n, 0)

        # Calculation of the first x, y curve positions
        x = c_nk * (self.t**0) * ((1 - self.t)**(n - 0)) * self.p[0][0]
        y = c_nk * (self.t**0) * ((1 - self.t)**(n - 0)) * self.p[0][1]
        z = c_nk * (self.t**0) * ((1 - self.t)**(n - 0)) * self.p[0][2]

        for i in range(1, n + 1):
            # Binomial cooficient in interation (i)
            c_nk = self.__binomial_coefficient(n, i)

            # Calculation positions in iteration (i)
            x_aux, y_aux, z_aux = self.__n_index_curve(i, self.p[i], self.t, n, c_nk)

            # The sum of all positions for the resulting Bézier curve
            x += x_aux
            y += y_aux
            z += z_aux

        self.__axis.plot3D(x, y, z, 'b--', label=r'N-Degree Bezier Curve: [$p_{0}$, $p_{1}$]; [$p_{1}$, $p_{2}$]; ... ; [$p_{n - 1}$, $p_{n}$]', linewidth=2.5)

    def __display_aux_result(self):
        """
        Description:
            Function for displaying points and text labeling.
        """

        for i in range(len(self.p)):
            self.__axis.text(self.p[i][0] + 0.01, self.p[i][1] + 0.01, self.p[i][2] + 0.01, '$p_{' + str(i) + '}$ = [' + str(self.p[i][0]) + ', ' + str(self.p[i][1]) + ', ' + str(self.p[i][2]) + ']', fontsize=20)

            if i != len(self.p) - 1:
                self.__axis.plot3D(self.p[i][0], self.p[i][1], self.p[i][2], marker = 'o', ms = 15, mfc = [1,1,1], markeredgecolor = [0,0,0], mew = 5)
            else:
                self.__axis.plot3D(self.p[i][0], self.p[i][1], self.p[i][2], label=self.__switch_dCtrl(len(self.p)), marker = 'o', ms = 15, mfc = [1,1,1], markeredgecolor = [0,0,0], mew = 5)

    def display_result(self):
        """
        Description:
            Function for calculating and displaying the results of Bézier curves.
        """

        try:
            assert len(self.p) > 1

            self.__display_aux_result()

            # If the number of user input points is greater than 4, the program chooses the n_points calculation method, but if the simplification 
            # coefficient is greater than 1, the program can choose another method and this calculation principle is faster.
            if self.simplification_factor > 1:
                # If the coefficient coefficient is greater than 1, simplify the path
                self.p = self.__path_simplification(self.p, self.simplification_factor)

            # Select a calculation method based on the number of points in the array (p).
            if len(self.p) == 2:
                self.__two_points()
            elif len(self.p) == 3:
                self.__three_points()
            elif len(self.p) == 4:
                self.__four_points()
            else:
                self.__n_points()

            # Set additional features for successful display of the Bézier curves.
            self.__plt.grid()
            self.__axis.set_xlabel('x axis [Unit]', fontsize = 20, fontweight ='normal')
            self.__axis.set_ylabel('y axis [Unit]', fontsize = 20, fontweight ='normal')
            self.__axis.set_zlabel('z axis [Unit]', fontsize = 20, fontweight ='normal')
            self.__plt.title('Bezier Curve', fontsize = 50, fontweight ='normal')
            self.__plt.legend(loc=0, fontsize=20)

            # Display a figure. Wait for the user to close the window.
            self.__plt.show()

        except AssertionError as error:
            print('[INFO] Insufficient number of entry points.')
            print('[INFO] The minimum number of entry points is 2.')


def main():
    # Initialization of the Class (Control Manipulator)
    # Input:
    #   (1 - 4) Points [Float Array]
    #   (5) Time Step  [INT]
    # Example:
    #   x = bezier_ctrl([1.0, 1.0], [1.25, 2.0], None, None, 100)

    # Try the calculation of the Bézier curve:
    # Select one of these options: 2 - Linear Curve, 3 - Quadratic Curve, 4 - Cubic Curve
    test = 'n_degree'

    if test == 'linear':
        # Linear Curve
        bezier = bezier_ctrl([[1.0, 1.0, 1.0], [1.25, 2.0, 2.0]], 100)
    elif test == 'quadratic':
        # Quadratic Curve
        bezier = bezier_ctrl([[1.0, 1.0, 1.0], [1.25, 2.0, 2.0], [1.75, 2.0, 1.5]], 100)
    elif test == 'cubic':
        # Cubic Curve
        bezier = bezier_ctrl([[1.0, 1.0, 1.0], [1.25, 2.0, 2.0], [1.75, 2.0, 1.5], [2.0, 1.0, 2.0]], 100)
    elif test == 'n_degree':
        bezier = bezier_ctrl([[1.0, 1.0, 1.0], [1.25, 2.0, 2.5], [1.75, 2.0, 1.5], [2.0, 1.0, 1.0], [1.0, -1.0, 2.0], [1.25, -2.0, 1.75], [1.75, -2.0, 2.75], [2.0, -1.0, 2.0]], 100)
        bezier.simplification_factor = 1

    # Display the result of the calculation -> figure with the resulting Bézier curves
    bezier.display_result()

if __name__ == '__main__':
    sys.exit(main())
