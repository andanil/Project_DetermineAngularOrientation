using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class Eigendecomp
    {
        public List<double> Eigenvalues { get; private set; }
        public List<Vector> Eigenvectors { get; private set; }

        public Eigendecomp()
        {
            Eigenvalues = new List<double>();
            Eigenvectors = new List<Vector>();
        }

        public void Compute(Matrix matrix, double eps)
        {
            Eigenvalues = ComputeEigenvalues(matrix, eps);
            Eigenvectors = new List<Vector>();
        }

        private List<double> ComputeEigenvalues(Matrix matrix, double eps)
        {
            double a = -1;
            double b = matrix[0, 0] + matrix[1, 1] + matrix[2, 2];
            double c = matrix[1, 0] * matrix[0, 1] + matrix[0, 2] * matrix[2, 0] + matrix[2, 1] * matrix[1, 2] -
                 matrix[0, 0] * matrix[1, 1] - matrix[0, 0] * matrix[2, 2] - matrix[1, 1] * matrix[2, 2];
            double d = matrix[0, 0] * matrix[1, 1] * matrix[2, 2] + matrix[1, 0] * matrix[2, 1] * matrix[0, 2] +
                matrix[2, 0] * matrix[0, 1] * matrix[1, 2] - matrix[0, 0] * matrix[1, 2] * matrix[2, 1] -
                matrix[1, 0] * matrix[0, 1] * matrix[2, 2] - matrix[2, 0] * matrix[1, 1] * matrix[0, 2];

            CubicEquation equation = new CubicEquation(a, b, c, d);
            List<double> eigenvalues = ToReal(equation.CardanosMethod(eps));
            eigenvalues.Sort();
            return eigenvalues;
        }

        private List<double> ToReal(List<Complex> list)
        {
            List<double> result = new List<double>();
            for (int i = 0; i < list.Count; i++)
                result.Add(list[i].Real);
            return result;
        }
    }
}
