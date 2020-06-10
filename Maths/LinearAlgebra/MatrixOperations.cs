using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class MatrixOperations
    {
        public static Eigendecomp Eigendecomposition(Matrix matrix, double eps)
        {
            Eigendecomp result = new Eigendecomp();
            result.Compute(matrix, eps);
            return result;
        }

        public static SVD SVD(Matrix matrix, double eps)
        {
            SVD result = new SVD();
            result.Compute(matrix, eps);
            return result;
        }

        public static Matrix RotationMatrix(double phi, double theta, double psi)
        {
            double[,] values = new double[,] {{ Math.Cos(phi)* Math.Cos(psi)- Math.Sin(phi)* Math.Sin(psi)* Math.Cos(theta),
                Math.Cos(phi) * Math.Sin(psi) + Math.Sin(phi) * Math.Cos(psi) * Math.Cos(theta), Math.Sin(phi) * Math.Sin(theta) },
            { -Math.Sin(phi)* Math.Cos(psi)- Math.Cos(phi)* Math.Sin(psi)* Math.Cos(theta),
                -Math.Sin(phi) * Math.Sin(psi) + Math.Cos(phi) * Math.Cos(psi) * Math.Cos(theta), Math.Cos(phi) * Math.Sin(theta) },
            { Math.Sin(theta) * Math.Sin(psi), -Math.Sin(theta) * Math.Cos(psi), Math.Cos(theta) } };

            return new Matrix(values);
        }
    }
}
