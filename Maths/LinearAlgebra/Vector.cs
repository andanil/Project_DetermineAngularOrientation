using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class Vector
    {
        private double[] values;
        private bool transpose;
        private int n = 3;

        public Vector()
        {
            values = new double[n];
            transpose = false;
        }

        public Vector( double[] value)
        {
            value.CopyTo(values, 0);
        }

        public double this[int i]
        {

            get { return values[i]; }
            set { values[i] = value; }
        }

        public double Nrm(int rowStart)
        {
            double s = 0;
            for (var i = rowStart; i < n; i++)
            {
                s += this[i] * this[i];
            }

            return Math.Sqrt(s);
        }

        public void Transpose()
        {
            transpose = !transpose;
        }

        public static Vector ZeroVector()
        {
            Vector vec = new Vector();
            vec.values = Enumerable.Repeat<double>(0, vec.n).ToArray();
            return vec;
        }
    }
}
