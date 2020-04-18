using System;
using Maths;
using Maths.LinearAlgebra;

namespace SvdSimpleApp
{
    class Program
    {
        static void Main(string[] args)
        {
            double[,] a1 = new double[,] { {0.53637763463723365, 0.028471325537409321, 0.66987587449600727 },
              { 0.96396422198226872, 0.82154775123183976, 0.84203247765173783 },
              { 0.30928549371160824, 0.36416072042852676, 0.19277459904215047 } };
           // Matrix a = Matrix.RandomMatrix3();
            SVD svd = MatrixOperations.SVD(new Matrix(a1), 1e-14);

            LibraryMethod(a1);

            Console.ReadKey();
        }

        static void EIG()
        {
            Matrix m1 = new Matrix(new double[,] { { 12, -34, 51 }, { -34, 40, 31 }, { 51, 31, -10 } });//три разных
            Matrix m2 = new Matrix(new double[,] { { 4, -5, 2 }, { 5, -7, 3 }, { 6, -9, 4 } }); //два одинаковых
            Matrix m3 = new Matrix(new double[,] { { 1, -3, 4 }, { 4, -7, 8 }, { 6, -7, 7 } }); //два одинаковых
            Matrix m4 = new Matrix(new double[,] { { 0, 1, 0 }, { -4, 4, 0 }, { -2, 1, 2 } }); //три равных
            Matrix m5 = new Matrix(new double[,] { { 1, 2, 1 }, { 2, -5, 2 }, { 1, 2, 1 } }); //три разных собственных значения
            Eigendecomp eigendecomp = MatrixOperations.Eigendecomposition(m2, 1e-14);
        }

        static void LibraryMethod(double[,] a)
        {
            double[,] u = new double[3, 3];
            double[,] vt = new double[3, 3];
            double[] w = new double[3];
            alglib.rmatrixsvd(a, 3, 3, 2, 2, 2, out w, out u, out vt);
            Console.WriteLine("-----------Library--------");
            Console.WriteLine("Matrix U");
            PrintData(u, 3, 3);
            Console.WriteLine("Matrix VT");
            PrintData(vt, 3, 3);
            Console.WriteLine("Vetor W");
            PrintData(w);
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
