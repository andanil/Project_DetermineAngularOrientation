using System;
using System.Collections.Generic;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class EigenDecomp
    {
        public List<double> EigenVal { get; private set; }
        public List<Vector> EigenVect { get; private set; }

        public EigenDecomp(List<double> values, List<Vector> vectors)
        {
            EigenVal = values;
            EigenVect = vectors;
        }

        public static EigenDecomp Compute(Matrix matrix)
        {
            List<double> eigenValues = new List<double>();
            List<Vector> eigenVectors = new List<Vector>();

            return new EigenDecomp(eigenValues, eigenVectors);
        }
    }
}
