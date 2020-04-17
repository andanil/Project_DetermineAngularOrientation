using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class EigenDecomp
    {
        public List<Complex> Eigenval { get; private set; }
        public List<Vector> Eigenvect { get; private set; }

        public EigenDecomp(List<Complex> values, List<Vector> vectors)
        {
            Eigenval = values;
            Eigenvect = vectors;
        }

        public static EigenDecomp Compute(Matrix matrix, double eps)
        {
            List<Complex> eigenvalues = ComputeEigenvalues(matrix, eps);
            List<Vector> eigenvectors = new List<Vector>();

            return new EigenDecomp(eigenvalues, eigenvectors);
        }

        private static List<Complex> ComputeEigenvalues(Matrix matrix, double eps)
        {
            double a = -1;
            double b = matrix[0,0] + matrix[1,1] + matrix[2,2];
            double c = matrix[1,0] * matrix[0,1] + matrix[0,2] * matrix[2,0] + matrix[2,1] * matrix[1,2] -
                 matrix[0, 0] * matrix[1, 1] - matrix[0, 0] * matrix[2, 2] - matrix[1, 1] * matrix[2, 2];
            double d = matrix[0, 0] * matrix[1, 1] * matrix[2, 2] + matrix[1, 0] * matrix[2, 1] * matrix[0, 2] +
                matrix[2, 0] * matrix[0, 1] * matrix[1, 2] - matrix[0, 0] * matrix[1, 2] * matrix[2, 1] -
                matrix[1, 0] * matrix[0, 1] * matrix[2, 2] - matrix[2, 0] * matrix[1, 1] * matrix[0, 2];

            CubicEquation equation = new CubicEquation(a, b, c, d);
            List<Complex> eigenvalues = equation.CardanosMethod(eps);
            eigenvalues.Sort();
            return eigenvalues;
        }

    }
}
