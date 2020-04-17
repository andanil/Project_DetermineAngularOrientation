using System;
using System.Collections.Generic;
using System.Linq;

namespace Maths.LinearAlgebra
{
    public class Matrix
    {
        private double[,] values;
        private int n = 3;

        public Matrix()
        {
            values = new double[n, n];
        }

        public Matrix(double [,] value)
        {
            value.CopyTo(values, 0);
        }

        public double this[int i, int j]
        {
            get { return values[i, j]; }
            set { values[i, j] = value; }
        }

        public Matrix Copy()
        {
            Matrix copy = new Matrix();
            for (int i = 0; i < n; i++)
                for(int j = 0; j < n; j++)
                copy.values[i, j] = values[i, j];
            return copy;
        }


        public static Matrix RandomMatrix3()
        {
            Matrix result = new Matrix();
            Random rnd = new Random();
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    result.values[i, j] = rnd.NextDouble();
            return result;
        }

        public Matrix Transpose()
        {
            Matrix trans = new Matrix();
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    trans[j, i] = this[i, j];
                }
            }
            return trans;
        }

        public static Matrix ZeroMatrix()
        {
            Matrix matrix = new Matrix();
            for (int i = 0; i < matrix.n; i++)
                for (int j = 0; j < matrix.n; j++)
                    matrix[i, j] = 0;
            return matrix;
        }

        public static Matrix IdentityMatrix()
        {
            Matrix matrix = new Matrix();
            for (int i = 0; i < matrix.n; i++)
                for (int j = 0; j < matrix.n; j++)
                    if(i == j)
                        matrix[i, j] = 1;
                    else
                        matrix[i, j] = 0;
            return matrix;
        }

        public static Matrix operator *(Matrix m1, Matrix m2)
        {
            Matrix r = Matrix.ZeroMatrix();
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        r[i, j] += m1[i, k] * m2[k, j];
            return r;
        }

        public static Vector operator *(Matrix m, Vector v)
        {
            Vector r = Vector.ZeroVector();
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    r[i] += m[i, j] * v[j];
            return r;
        }
        public static Vector operator *(Vector v, Matrix m)
        {
            Vector r = Vector.ZeroVector();
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    r[i] += m[j, i] * v[j];
            r.Transpose();
            return r;
        }

        public static Matrix operator *(double constant, Matrix m)
        {
            Matrix r = m.Copy();
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    r[i,j] = constant * m[i,j];
            return r;
        }
    }
}
