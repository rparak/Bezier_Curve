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
File Name: test.py
## =========================================================================== ## 
"""

# System (Default Lib.)
import sys

# Numpy (Array computing Lib.) [pip3 install numpy]
import numpy as np

# Mtaplotlib (Visualization Lib.) [pip3 install matplotlib]
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

# Own library for Bézier curve calculation
import Bezier

def main():
    # Select one of these options: '2D', '3D' {Dimensional Space}
    dimensional_space = '3D'
    
    # Visibility of curves
    visible_linear, visible_quadratic, visible_cubic, visible_auto = True, True, True, True
    
    if dimensional_space == '2D':
        points = [[1.0, 1.0], [1.25, 2.0], [1.75, 2.0], [2.0, 1.0], [1.0, -1.0], [1.25, -2.0], [1.75, -2.0], [2.0, -1.0]]
    else:
        _, axis = plt.subplots(subplot_kw={"projection": "3d"})
        points = [[1.0, 1.0, 1.0], [1.25, 2.0, 2.5], [1.75, 2.0, 1.5], [2.0, 1.0, 1.0], [1.0, -1.0, 2.0], [1.25, -2.0, 1.75], [1.75, -2.0, 2.75], [2.0, -1.0, 2.0]]

    """
    Initialization Bezier Class
    
    Input:
        (1) Time Step  [INT]
    """
    ts          = 100
    Bezier_Ndeg = Bezier.N_Degree(ts)

    # Time t ∈ [0, 1]
    time = np.linspace(0.0, 1.0, ts)

    """
    Description:
        1\ Calculation and display the Linear Bézier Curve p(t).
    """
    if visible_linear == True:
        for i in range(len(points) - 1):
            result = Bezier.Linear(time, [points[i], points[i + 1]])
            
            if i == (len(points) - 1) - 1:
                if dimensional_space == '2D':
                    plt.plot(result[0], result[1], 'r--', label=r'Linear Bezier Curve: [$p_{0}$, $p_{1}$] ; ... ; [$p_{n - 1}$, $p_{n}$]', linewidth=2.5)
                else:
                    axis.plot(result[0], result[1], result[2], 'r--', label=r'Linear Bezier Curve: [$p_{0}$, $p_{1}$] ; ... ; [$p_{n - 1}$, $p_{n}$]' , linewidth=2.5)
            else:
                if dimensional_space == '2D':
                    plt.plot(result[0], result[1], 'r--', linewidth=2.5)
                else:
                    axis.plot(result[0], result[1], result[2], 'r--', linewidth=2.5)

    """
    Description:
        2\ Calculation and display the Quadratic Bézier Curve p(t).
    """
    if visible_quadratic == True:
        for i in range(len(points) - 2):
            result = Bezier.Quadratic(time, [points[i], points[i + 1], points[i + 2]])
            
            if i == (len(points) - 2) - 1:
                if dimensional_space == '2D':
                    plt.plot(result[0], result[1], 'g--', label=r'Quadratic Bezier Curve: [$p_{0}$, $p_{1}$, $p_{2}$] ; ... ; [$p_{n - 2}$, $p_{n - 1}$, $p_{n}$]', linewidth=2.5)
                else:
                    axis.plot(result[0], result[1], result[2], 'g--', label=r'Quadratic Bezier Curve: [$p_{0}$, $p_{1}$, $p_{2}$] ; ... ; [$p_{n - 2}$, $p_{n - 1}$, $p_{n}$]' , linewidth=2.5)
            else:
                if dimensional_space == '2D':
                    plt.plot(result[0], result[1], 'g--', linewidth=2.5)
                else:
                    axis.plot(result[0], result[1], result[2], 'g--', linewidth=2.5)

    """
    Description:
        3\ Calculation and display the Cubic Bézier Curve p(t).
    """
    if visible_cubic == True:
        for i in range(len(points) - 3):
            result = Bezier.Cubic(time, [points[i], points[i + 1], points[i + 2], points[i + 3]])
            
            if i == (len(points) - 3) - 1:
                if dimensional_space == '2D':
                    plt.plot(result[0], result[1], 'b--', label=r'Cubic Bezier Curve: [$p_{0}$, $p_{1}$, $p_{2}$, $p_{3}$] ; ... ; [$p_{n - 3}$, $p_{n - 2}$, $p_{n - 1}$, $p_{n}$]', linewidth=2.5)
                else:
                    axis.plot(result[0], result[1], result[2], 'b--', label=r'Cubic Bezier Curve: [$p_{0}$, $p_{1}$, $p_{2}$, $p_{3}$] ; ... ; [$p_{n - 3}$, $p_{n - 2}$, $p_{n - 1}$, $p_{n}$]' , linewidth=2.5)
            else:
                if dimensional_space == '2D':
                    plt.plot(result[0], result[1], 'b--', linewidth=2.5)
                else:
                    axis.plot(result[0], result[1], result[2], 'b--', linewidth=2.5)

    """
    Description:
        4\ Automatic Curve Calculation (Result of the calculation) -> depends on the simplification factor
    """
    if visible_auto == True:
        result = Bezier_Ndeg.Solve(points, 1)

        # Display the N-Degree Bézier Curve p(t)
        if dimensional_space == '2D':
            plt.plot(result[0], result[1], 'm--', label=r'N-Degree Bezier Curve: [$p_{0}$, $p_{1}$, ... , $p_{n - 1}$, $p_{n}$]', linewidth=2.5)
        elif dimensional_space == '3D':
            axis.plot(result[0], result[1], result[2], 'm--', label=r'N-Degree Bezier Curve: [$p_{0}$, $p_{1}$, ... , $p_{n - 1}$, $p_{n}$]' , linewidth=2.5)

    if dimensional_space == '2D':
        for _, point in enumerate(points):
            plt.plot(point[0], point[1], marker = 'o', ms = 15, mfc = [1,1,1], markeredgecolor = [0,0,0], mew = 5)
    elif dimensional_space == '3D':
        for _, point in enumerate(points):
            axis.scatter(point[0], point[1], point[2], marker = 'o', color = [1,1,1], edgecolors = [0,0,0], linewidths=2.5)

    # Set additional features for successful display of the Bézier curves.
    plt.grid()
    plt.xlabel('x axis [Unit]', fontsize = 20, fontweight ='normal')
    plt.ylabel('y axis [Unit]', fontsize = 20, fontweight ='normal')
    plt.title('Bezier Curve', fontsize = 50, fontweight ='normal')
    plt.legend(loc=0,fontsize=20)

    # Display the result of the calculation -> figure with the resulting Bézier curve
    plt.show()

if __name__ == '__main__':
    sys.exit(main())