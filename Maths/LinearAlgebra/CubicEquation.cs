using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class CubicEquation
    {
        private double a;
        private double b;
        private double c;
        private double d;
        public CubicEquation(double a, double b, double c, double d)
        {
            this.a = a;
            this.b = b;
            this.c = c;
            this.d = d;
        }

        public  List<Complex> CardanosMethod(double eps)
        {
            double alpha;
            double betta;
            List<Complex> x = new List<Complex>();
            Normalize();
            double p = - b * b / 3 + c;
            double q = 2 * Math.Pow(b, 3) / 27 - b * c / 3 *  + d;
            double D = q * q / 4 + Math.Pow(p, 3) / 27;
            if(Math.Abs(D) < eps)
            {
                alpha = CubeRoot(-q / 2);
                x.Add(new Complex(2 * alpha, 0));
                x.Add(new Complex(-alpha, 0));
                x.Add(new Complex(-alpha, 0));
            }
            else
            {
                if(D > 0)
                {
                    alpha = CubeRoot(-q / 2 + Math.Sqrt(D));
                    betta = CubeRoot(-q / 2 - Math.Sqrt(D));
                    x[0] = new Complex(alpha + betta, 0);
                    x[1] = new Complex(-(alpha + betta) / 2, Math.Sqrt(3)*(alpha - betta)/2);
                    x[2] = new Complex(x[1].Real, -x[1].Imaginary);
                }
                else
                {
                    Complex buf = new Complex(-q / 2, Math.Sqrt(-D));
                    double r = buf.Magnitude;
                    r = CubeRoot(r);
                    double phi = buf.Phase;
                    x[0] = new Complex(2 * r * Math.Cos(phi / 3), 0);
                    x[1] = new Complex(2 * r * Math.Cos((phi + 2 * Math.PI) / 3), 0);
                    x[2] = new Complex(2 * r * Math.Cos((phi + 4 * Math.PI) / 3), 0);
                }
            }
            for (int i = 0; i < 3; i++)
                x[i] -= b / 3;
            return x;
        }
        
        private double CubeRoot(double x)
        {
            if (x >= 0)
                return Math.Pow(x, (1 / 3));
            return -Math.Pow(-x, (1 / 3));
        }

        private void Normalize()
        {
            b /= a;
            c /= a;
            d /= a;
        }
    }
}
