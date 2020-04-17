using System;
using Maths;
using Maths.LinearAlgebra;

namespace SvdSimpleApp
{
    class Program
    {
        static void Main(string[] args)
        {
            /*double[,] a = new double[,] { {0.29155327113883256, 0.74628977000074914, 0.29136221403785151 },
              { 0.928606090102627, 0.6366486431270133, 0.223203904565053 },
              { 0.84984116808038257, 0.57135150841034554, 0.91516676541192776 } };
            double[,] u = new double[3, 3];
            double[,] vt = new double[3, 3];
            double[] w = new double[3];
            alglib.rmatrixsvd(a, 3, 3, 2, 2, 2, out w, out u, out vt);
            Console.WriteLine("Matrix U");
            PrintData(u, 3, 3);
            Console.WriteLine("Matrix VT");
            PrintData(vt, 3, 3);
            Console.WriteLine("Vetor W");
            PrintData(w);*/
            Eigendecomp eigendecomp = MatrixOperations.Eigendecomposition(new Matrix(new double[,] { { 12, -34, 51 }, { -34, 40, 31 }, { 51, 31, -10 } }), 1e-14);

            Console.ReadKey();
        }

        static void PrintData(double[,] values, int n, int m)
        {
            for(int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                    Console.Write(values[i, j].ToString() + " ");
                Console.Write("\n");
            }
        }

        static void PrintData(double[] values)
        {
            for (int i = 0; i < values.Length; i++)
                Console.Write(values[i].ToString() + " ");
        }
    }
}
