using System;
using System.Collections.Generic;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class SVD
    {
        public Matrix S { get; protected set; }
        public Matrix V { get; protected set; }
        public Matrix U { get; protected set; }

        public SVD()
        {
            S = new Matrix();
            V = new Matrix();
            U = new Matrix();
        }

        public void Compute(Matrix matrix, double eps)
        {
            Matrix g = matrix.Transpose() * matrix;
            Eigendecomp eigendecomp = MatrixOperations.Eigendecomposition(g, eps);
            SortEigen(eigendecomp);
            Vector sigma = new Vector();
            for (int i = 0; i < matrix.N; i++)
                sigma[i] = Math.Sqrt(Math.Abs(eigendecomp.Eigenvalues[i]));
            V = new Matrix();
            S = Matrix.ZeroMatrix();
            U = new Matrix(eigendecomp.Eigenvectors);
            List<Vector> v = new List<Vector>();
            for(int i = 0; i < matrix.N; i++)
            {
                S[i, i] = sigma[i];
                v.Add(matrix * U.GetColumn(i)/sigma[i]);
            }
            V = new Matrix(v);

        }

        public double Error(Matrix matrix)
        {
            Matrix matrix1 = V * S * U.Transpose();
            Matrix result = matrix - matrix1;
            return result.InfinityNorm();
        }

        private void SortEigen(Eigendecomp eigendecomp)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3 - i - 1; j++)
                    if (eigendecomp.Eigenvalues[j] < eigendecomp.Eigenvalues[j + 1])
                        Swap(eigendecomp, j, j + 1);
        }

        private void Swap(Eigendecomp eigendecomp, int first, int second)
        {
            double b = eigendecomp.Eigenvalues[first];
            eigendecomp.Eigenvalues[first] = eigendecomp.Eigenvalues[second];
            eigendecomp.Eigenvalues[second] = b;
            Vector vec = eigendecomp.Eigenvectors[first].Copy();
            eigendecomp.Eigenvectors[first] = eigendecomp.Eigenvectors[second].Copy();
            eigendecomp.Eigenvectors[second] = vec;
        }
    }
}
