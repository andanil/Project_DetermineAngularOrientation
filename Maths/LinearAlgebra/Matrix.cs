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
            values = value;
        }

        public Matrix(List<Vector> vectors)
        {
            values = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    values[i, j] = vectors[j][i];
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
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    trans[j, i] = this[i, j];
                }
            }
            return trans;
        }

        public Vector GetRow(int index)
        {
            Vector row = new Vector();
            for (int i = 0; i < n; i++)
                row[i] = this[index, i];
            return row;
        }

        public Vector GetColumn(int index)
        {
            Vector column = new Vector();
            for (int i = 0; i < n; i++)
                column[i] = this[i, index];
            return column;
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
            Matrix result = ZeroMatrix();
            for (int i = 0; i < m1.n; i++)
                for (int j = 0; j < m1.n; j++)
                    for (int k = 0; k < m1.n; k++)
                        result[i, j] += m1[i, k] * m2[k, j];
            return result;
        }

        public static Vector operator *(Matrix m, Vector v)
        {
            Vector result = Vector.ZeroVector();
            for (int i = 0; i < m.n; i++)
                for (int j = 0; j < m.n; j++)
                    result[i] += m[i, j] * v[j];
            return result;
        }
        public static Vector operator *(Vector v, Matrix m)
        {
            Vector result = Vector.ZeroVector();
            for (int i = 0; i < m.n; i++)
                for (int j = 0; j < m.n; j++)
                    result[i] += m[j, i] * v[j];
            result.Transpose();
            return result;
        }

        public static Matrix operator *(double constant, Matrix m)
        {
            Matrix result = m.Copy();
            for (int i = 0; i < m.n; i++)
                for (int j = 0; j < m.n; j++)
                    result[i,j] = constant * m[i,j];
            return result;
        }
    }
}
