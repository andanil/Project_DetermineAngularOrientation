using System;
using System.Collections.Generic;
using System.Linq;

namespace Maths.LinearAlgebra
{
    public class Matrix
    {
        private double[,] values;

        public int N { get; private set; } = 3;

        public Matrix()
        {
            values = new double[N, N];
        }

        public Matrix(double [,] value)
        {
            values = value;
        }

        public Matrix(List<Vector> vectors)
        {
            values = new double[N, N];
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    values[i, j] = vectors[j][i];
        }

        public double this[int i, int j]
        {
            get { return values[i, j]; }
            set { values[i, j] = value; }
        }

        public double InfinityNorm()
        {
            double norm = 0;
            double sum;
            for(int i = 0; i < N; i++)
            {
                sum = 0;
                for (int j = 0; j < N; j++)
                    sum += Math.Abs(values[i, j]);
                if (sum > norm)
                    norm = sum;
            }
            return norm;
        }

        public Matrix Copy()
        {
            Matrix copy = new Matrix();
            for (int i = 0; i < N; i++)
                for(int j = 0; j < N; j++)
                copy.values[i, j] = values[i, j];
            return copy;
        }

        public double[,] ToArray()
        {
            return values;
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
            for (int j = 0; j < N; j++)
            {
                for (int i = 0; i < N; i++)
                {
                    trans[j, i] = values[i, j];
                }
            }
            return trans;
        }

        public Vector GetRow(int index)
        {
            Vector row = new Vector();
            for (int i = 0; i < N; i++)
                row[i] = values[index, i];
            return row;
        }

        public Vector GetColumn(int index)
        {
            Vector column = new Vector();
            for (int i = 0; i < N; i++)
                column[i] = values[i, index];
            return column;
        }

        public static Matrix ZeroMatrix()
        {
            Matrix matrix = new Matrix();
            for (int i = 0; i < matrix.N; i++)
                for (int j = 0; j < matrix.N; j++)
                    matrix[i, j] = 0;
            return matrix;
        }

        public static Matrix IdentityMatrix()
        {
            Matrix matrix = new Matrix();
            for (int i = 0; i < matrix.N; i++)
                for (int j = 0; j < matrix.N; j++)
                    if(i == j)
                        matrix[i, j] = 1;
                    else
                        matrix[i, j] = 0;
            return matrix;
        }

        public static Matrix operator *(Matrix m1, Matrix m2)
        {
            Matrix result = ZeroMatrix();
            for (int i = 0; i < m1.N; i++)
                for (int j = 0; j < m1.N; j++)
                    for (int k = 0; k < m1.N; k++)
                        result[i, j] += m1[i, k] * m2[k, j];
            return result;
        }

        public static Vector operator *(Matrix m, Vector v)
        {
            Vector result = Vector.ZeroVector();
            for (int i = 0; i < m.N; i++)
                for (int j = 0; j < m.N; j++)
                    result[i] += m[i, j] * v[j];
            return result;
        }
        public static Vector operator *(Vector v, Matrix m)
        {
            Vector result = Vector.ZeroVector();
            for (int i = 0; i < m.N; i++)
                for (int j = 0; j < m.N; j++)
                    result[i] += m[j, i] * v[j];
            result.Transpose();
            return result;
        }

        public static Matrix operator *(double constant, Matrix m)
        {
            Matrix result = m.Copy();
            for (int i = 0; i < m.N; i++)
                for (int j = 0; j < m.N; j++)
                    result[i,j] = constant * m[i,j];
            return result;
        }

        public static Matrix operator -(Matrix m1, Matrix m2)
        {
            Matrix result = ZeroMatrix();
            for (int i = 0; i < m1.N; i++)
                for (int j = 0; j < m1.N; j++)
                        result[i, j] = m1[i, j] - m2[i, j];
            return result;
        }
    }
}
