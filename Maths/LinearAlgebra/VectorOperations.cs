using System;
using System.Collections.Generic;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class VectorOperations
    {
        public static Vector Cross(Vector v1, Vector v2)
        {
            Vector vector = new Vector();
            vector[0] = v1[1] * v2[2] - v1[2] * v2[1];
            vector[1] = v1[2] * v2[0] - v1[0] * v2[2];
            vector[2] = v1[0] * v2[1] - v1[1] * v2[0];
            return vector;
        }

        public static void OrthogonalComplement(Vector w, out Vector v, out Vector u)
        {
            double invLength;
            if(Math.Abs(w[0]) > Math.Abs(w[1]))
            {
                invLength = 1 / Math.Sqrt(w[0] * w[0] + w[2] * w[2]);
                u = new Vector(new double[] { -w[2] * invLength, 0, w[0] * invLength });
            }
            else
            {
                invLength = 1 / Math.Sqrt(w[1] * w[1] + w[2] * w[2]);
                u = new Vector(new double[] { 0, w[2] * invLength, -w[1] * invLength });
            }
            v = Cross(u, w);
        }
    }
}
