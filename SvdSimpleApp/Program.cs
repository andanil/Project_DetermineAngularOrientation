using System;
using System.Linq;
using System.Diagnostics;
using Maths;
using Maths.LinearAlgebra;
using System.Collections.Generic;

namespace SvdSimpleApp
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.BackgroundColor = ConsoleColor.White;
            Console.ForegroundColor = ConsoleColor.Black;
            Console.Clear();
            /*double[,] a1 = new double[,] { {0.53637763463723365, 0.028471325537409321, 0.66987587449600727 },
              { 0.96396422198226872, 0.82154775123183976, 0.84203247765173783 },
              { 0.30928549371160824, 0.36416072042852676, 0.19277459904215047 } };*/
            Matrix a = Matrix.RandomMatrix3();
            // Matrix a = new Matrix(a1);
            Console.WriteLine("Матрица A");
            PrintData(a.ToArray(), 3, 3);

            SVD svd = MatrixOperations.SVD(a, 1e-15);
            Console.WriteLine();
            Console.WriteLine("Матрица V");
            PrintData(svd.V.ToArray(), 3, 3);
            Console.WriteLine();
            Console.WriteLine("Матрица S");
            PrintData(svd.S.ToArray(), 3, 3);
            Console.WriteLine();
            Console.WriteLine("Матрица U");
            PrintData(svd.U.ToArray(), 3, 3);

            double result = svd.Error(a);

            Console.WriteLine();
            Console.WriteLine("Погрешность:" + result.ToString());

            /*Stopwatch stopwatch = Stopwatch.StartNew();
            LibraryMethod(Random(1000, 1000), 1000, 1000);
            stopwatch.Stop();
            TimeSpan time = stopwatch.Elapsed;*/

            //Eigendecomp eigendecomp = EIG();

            /*List<Vector> vectors = new List<Vector>();
            vectors.Add(new Vector(new double[] { 0, -80, 0 }));
            vectors.Add(new Vector(new double[] { -69.3, 40, 0 }));
            vectors.Add(new Vector(new double[] { 69.3, 40, 0 }));
            vectors.Add(new Vector(new double[] { 0, 150, 0 }));
            vectors.Add(new Vector(new double[] { 0, 200, 0 }));

            List<Matrix> matrices = SimpleExperiments.Experiment(vectors, new Vector(new double[] { 50, -30, 5 }), 1e-13);
            Console.WriteLine("Матрица R");
            PrintData(matrices[0].ToArray(), 3, 3);
            Console.WriteLine("Матрица R'");
            PrintData(matrices[1].ToArray(), 3, 3);

            Console.WriteLine();
            Console.Write("Разница равномерных норм: ");
            Console.WriteLine(Math.Abs(matrices[0].InfinityNorm() - matrices[1].InfinityNorm()).ToString());*/

            Console.ReadKey();
        }

        static Eigendecomp EIG()
        {
            Matrix m1 = new Matrix(new double[,] { { 1, 2, 1 }, { 2, -5, 2 }, { 1, 2, 1 } });//три разных Проскуряков 1253
            Matrix m2 = new Matrix(new double[,] { { 1, 2, 2 }, { 2, 1, 2 }, { 2, 2, 1 } }); //два одинаковых Проскуряков 1251
            Matrix m3 = new Matrix(new double[,] { { 17, -2, -2 }, { -2, 14, -4 }, { -2, -4, 14 } }); //два одинаковых Проскуряков 1252
            Matrix m4 = new Matrix(new double[,] { { 8, 4, -1 }, { 4, -7, 4 }, { -1, 4, 8 } }); //два одинаковых Проскуряков 1254
            Matrix m5 = new Matrix(new double[,] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } }); 
            return MatrixOperations.Eigendecomposition(m5, 1e-14);
        }

        static void LibraryMethod(double[,] a, int n, int m)
        {
            double[,] u = new double[n, n];
            double[,] vt = new double[m, m];
            double[] w = new double[m];
            alglib.rmatrixsvd(a, n, m, 2, 2, 2, out w, out u, out vt);
            /*Console.WriteLine("-----------Library--------");
            Console.WriteLine("Matrix U");
            PrintData(u, n, 3);
            Console.WriteLine("Matrix VT");
            PrintData(vt, 3, 3);
            Console.WriteLine("Vetor W");
            PrintData(w);*/
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

        static double[,] Random(int n, int m)
        {
            double[,] result = new double[n,m];
            Random rnd = new Random();
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                    result[i, j] = rnd.NextDouble();
            return result;
        }
    }
}
