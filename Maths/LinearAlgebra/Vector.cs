using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class Vector
    {
        private double[] values;
        private int n = 3;

        public bool IsTransposed { get; private set; }

        public Vector()
        {
            values = new double[n];
            IsTransposed = false;
        }

        public Vector( double[] value)
        {
            values = value;
        }

        public Vector Copy()
        {
            return new Vector(values);
        }

        public double this[int i]
        {

            get { return values[i]; }
            set { values[i] = value; }
        }

        public double Norm()
        {
            double s = 0;
            for (var i = 0; i < n; i++)
            {
                s += this[i] * this[i];
            }
            return Math.Sqrt(s);
        }

        public void Transpose()
        {
            IsTransposed = !IsTransposed;
        }

        public static Vector ZeroVector()
        {
            Vector result = new Vector();
            result.values = Enumerable.Repeat<double>(0, result.n).ToArray();
            return result;
        }

        public static double operator *(Vector v1, Vector v2)
        {
            double result = 0;
            for (int i = 0; i < v1.n; i++)
                result += v1[i] * v2[i];
            return result;
        }

        public static Vector operator *(double constant, Vector v)
        {
            Vector result = new Vector();
            for (int i = 0; i < v.n; i++)
                result[i] = v[i] * constant;
            return result;
        }

        public static Vector operator *(Vector v, double constant)
        {
            Vector result = new Vector();
            for (int i = 0; i < v.n; i++)
                result[i] = v[i] * constant;
            return result;
        }

        public static Vector operator /(Vector v, double constant)
        {
            Vector result = new Vector();
            for (int i = 0; i < v.n; i++)
                result[i] = v[i] / constant;
            return result;
        }

        public static Vector operator -(Vector v1, Vector v2)
        {
            Vector result = new Vector();
            for (int i = 0; i < v1.n; i++)
                result[i] = v1[i] - v2[i];
            return result;
        }
    }
}
