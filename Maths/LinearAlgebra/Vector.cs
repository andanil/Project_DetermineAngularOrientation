using System;
using System.Collections.Generic;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class Vector
    {
        private double[] values;
        private bool transpose;

        public Vector()
        {
            values = new double[3];
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
            for (var i = rowStart; i < 3; i++)
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
            return new Vector(new double[] { 0, 0, 0 });
        }
    }
}
