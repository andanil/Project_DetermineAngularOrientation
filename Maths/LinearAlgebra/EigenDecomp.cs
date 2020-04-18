using System;
using System.Collections.Generic;
using System.Linq;
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
            Eigenvectors = ComputeEigenvectors(matrix, Eigenvalues, eps);
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

        private List<Vector> ComputeEigenvectors(Matrix matrix, List<double> eigenvalues, double eps)
        {
            List<Vector> eigenvectors = new List<Vector>() { new Vector(), new Vector(), new Vector() };
            if(Math.Abs(eigenvalues[0] - eigenvalues[1]) > eps)
            {
                eigenvectors[0] = ComputeFirstEigenvector(matrix, eigenvalues[0]);
                eigenvectors[1] = ComputeSecondEigenvector(matrix, eigenvalues[1], eigenvectors[0]);
                eigenvectors[2] = VectorOperations.Cross(eigenvectors[0], eigenvectors[1]);
            }
            else
            {
                eigenvectors[2] = ComputeFirstEigenvector(matrix, eigenvalues[2]);
                eigenvectors[1] = ComputeSecondEigenvector(matrix, eigenvalues[1], eigenvectors[2]);
                eigenvectors[0] = VectorOperations.Cross(eigenvectors[2], eigenvectors[1]);
            }

            /*for(int i = 0; i < eigenvalues.Count; i++)
                eigenvectors[i] = ComputeFirstEigenvector(matrix, eigenvalues[i]);*/

            return eigenvectors;
        }

        private Vector ComputeFirstEigenvector(Matrix matrix, double eigenvalue)
        {
            List<Vector> rows = new List<Vector>();
            for (int i = 0; i < 3; i++)
            {
                rows.Add(matrix.GetRow(i));
                rows[i][i] -= eigenvalue;
            }

            List<Vector> crossProd = new List<Vector>() { VectorOperations.Cross(rows[0], rows[1]),
                VectorOperations.Cross(rows[0], rows[2]), VectorOperations.Cross(rows[2], rows[1])};
            List<double> crossProdNrm = new List<double>();
            for (int i = 0; i < 3; i++)
                crossProdNrm.Add(crossProd[i].Nrm());
            int indexOfMax = crossProdNrm.IndexOf(crossProdNrm.Max());

            return crossProd[indexOfMax]/crossProdNrm[indexOfMax];
        }

        private Vector ComputeSecondEigenvector(Matrix matrix, double eigenvalue, Vector firstEigenvector)
        {
            Vector u, v, eigenvector;
            VectorOperations.OrthogonalComplement(firstEigenvector, out v, out u);

            Matrix a = matrix.Copy();
            for (int i = 0; i < 3; i++)
                a[i, i] -= eigenvalue;

            Vector ua = u * a;
            Vector va = v * a;
            double[,] m = new double[,] { { u * ua, v * ua }, { u * va, v * va} };
            double[,] absM = new double[2,2];
            for (int i = 0; i < 2; i++)
                for(int j = 0; j < 2; j++)
                    absM[i,j] = Math.Abs(m[i,j]);

            if(absM[0,0] >= absM[1,1])
            {
                double maxAbs = Math.Max(absM[0,0], absM[0,1]);
                if(maxAbs > 0)
                {
                    if(absM[0,0] >= absM[0,1])
                    {
                        m[0, 1] /= m[0, 0];
                        m[0, 0] = 1 / Math.Sqrt(1 + m[0, 1] * m[0, 1]);
                        m[0, 1] *= m[0, 0];
                    }
                    else
                    {
                        m[0, 0] /= m[0, 1];
                        m[0, 1] = 1 / Math.Sqrt(1 + m[0, 0] * m[0, 0]);
                        m[0, 0] *= m[0, 1];
                    }
                    eigenvector = m[0,1] * u - m[0,0] * v;
                }
                else
                {
                    eigenvector = u;
                }
            }
            else
            {
                double maxAbs = Math.Max(absM[1,1], absM[1, 0]);
                if (maxAbs > 0)
                {
                    if (absM[1, 1] >= absM[1,0])
                    {
                        m[1, 0] /= m[1, 1];
                        m[1, 1] = 1 / Math.Sqrt(1 + m[1, 0] * m[1, 0]);
                        m[1, 0] *= m[1, 1];
                    }
                    else
                    {
                        m[1, 1] /= m[1, 0];
                        m[1, 0] = 1 / Math.Sqrt(1 + m[1, 1] * m[1, 1]);
                        m[1, 1] *= m[1, 0];
                    }
                    eigenvector = m[1, 1] * u - m[1, 0] * v;
                }
                else
                {
                    eigenvector = u;
                }
            }
            
            return eigenvector;
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
