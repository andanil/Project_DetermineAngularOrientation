using System;
using System.Collections.Generic;
using System.Text;
using Maths.LinearAlgebra;

namespace SvdSimpleApp
{
    public class SimpleExperiments
    {
        public static List<Matrix> Experiment(List<Vector> vectors, Vector point, double eps)
        {
            List<Matrix> result = new List<Matrix>();
            Random rnd = new Random();
            double phi = rnd.NextDouble() * 360 - 180;
            double psi = rnd.NextDouble() * 60 - 30;
            double theta = rnd.NextDouble() * 60 - 30;
            Matrix rotMatrix = MatrixOperations.RotationMatrix(phi, theta, psi);
            result.Add(rotMatrix);

            Matrix E = ComputeEMatrix(vectors, point);
            Matrix moveE = rotMatrix * E;

            SVD svd = MatrixOperations.SVD(E*moveE.Transpose(), eps);
            Matrix R = svd.U * svd.V.Transpose();

            result.Add(R);
            return result;
        }

        private static Matrix ComputeEMatrix(List<Vector> vectors, Vector point)
        {
            List<Vector> e = new List<Vector>();
            for(int i = 0; i < vectors.Count; i++)
                e.Add((vectors[i] - point) / (vectors[i] - point).Norm());
            return new Matrix(e);
        }
    }
}
