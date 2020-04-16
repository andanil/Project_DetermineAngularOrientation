using System;
using System.Collections.Generic;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class SVD
    {
        public Vector S { get; protected set; }
        public Matrix V { get; protected set; }
        public Matrix U { get; protected set; }

        private SVD(Vector s, Matrix v, Matrix u)
        {
            S = s;
            V = v;
            U = u;
        }

        /*public Matrix ComputeW()
        {
            Matrix w = MatrixOperations.ZeroMatrix(U.RowCount, V.ColumnCount);
            for (int i = 0; i < S.Count; i++)
                w[i, i] = S[i];
            return w;
        }

        public static SVD Compute(Matrix matrix, int maxiter)
        {
            int nm = Math.Min(matrix.RowCount + 1, matrix.ColumnCount);
            Matrix matrixCopy = matrix.Copy();

            Vector s = new Vector(nm);
            Matrix u = new Matrix(matrixCopy.RowCount, matrixCopy.RowCount);
            Matrix vt = new Matrix(matrixCopy.ColumnCount, matrixCopy.ColumnCount);

            Vector e = new Vector(matrixCopy.ColumnCount);
            var work = new double[matrixCopy.RowCount];

            int i, j;
            int l, lp1;
            double t;

            var ncu = matrixCopy.RowCount;

            // Reduce matrixCopy to bidiagonal form, storing the diagonal elements
            // In s and the super-diagonal elements in e.
            var nct = Math.Min(matrixCopy.RowCount - 1, matrixCopy.ColumnCount);
            var nrt = Math.Max(0, Math.Min(matrixCopy.ColumnCount - 2, matrixCopy.RowCount));
            var lu = Math.Max(nct, nrt);
            for (l = 0; l < lu; l++)
            {
                lp1 = l + 1;
                if (l < nct)
                {
                    // Compute the transformation for the l-th column and place the l-th diagonal in VectorS[l].
                    double xnorm = NrmColumn(matrixCopy, matrixCopy.RowCount, l, l);
                    s[l] = xnorm;
                    if (s[l] != 0.0)
                    {
                        if (matrixCopy[l, l] != 0.0)
                        {
                            s[l] = Dsign(s[l], matrixCopy[l, l]);
                        }

                        ScalColumn(matrixCopy, matrixCopy.RowCount, l, l, 1.0 / s[l]);
                        matrixCopy[l, l] = 1.0 + matrixCopy[l, l];
                    }

                    s[l] = -s[l];
                }

                for (j = lp1; j < matrixCopy.ColumnCount; j++)
                {
                    if (l < nct)
                    {
                        if (s[l] != 0.0)
                        {
                            // Apply the transformation.
                            t = -DotColumns(matrixCopy, matrixCopy.RowCount, l, j, l) / matrixCopy[l, l];
                            for (var ii = l; ii < matrixCopy.RowCount; ii++)
                            {
                                matrixCopy[ii, j] = matrixCopy[ii, j] + t * matrixCopy[ii, l];
                            }
                        }
                    }

                    // Place the l-th row of matrixCopy into  e for the
                    // Subsequent calculation of the row transformation.
                    e[j] = matrixCopy[l, j];
                }

                if (l < nct)
                {
                    // Place the transformation in u for subsequent back multiplication.
                    for (i = l; i < matrixCopy.RowCount; i++)
                    {
                        u[i, l] = matrixCopy[i, l];
                    }
                }

                if (l >= nrt)
                {
                    continue;
                }

                // Compute the l-th row transformation and place the l-th super-diagonal in e(l).
                var enorm = e.Nrm(lp1);
                e[l] = enorm;
                if (e[l] != 0.0)
                {
                    if (e[lp1] != 0.0)
                    {
                        e[l] = Dsign(e[l], e[lp1]);
                    }

                    e.MultScal(lp1, 1.0 / e[l]);
                    e[lp1] = 1.0 + e[lp1];
                }

                e[l] = -e[l];
                if (lp1 < matrixCopy.RowCount && e[l] != 0.0)
                {
                    // Apply the transformation.
                    for (i = lp1; i < matrixCopy.RowCount; i++)
                    {
                        work[i] = 0.0;
                    }

                    for (j = lp1; j < matrixCopy.ColumnCount; j++)
                    {
                        for (var ii = lp1; ii < matrixCopy.RowCount; ii++)
                        {
                            work[ii] += e[j] * matrixCopy[ii, j];
                        }
                    }

                    for (j = lp1; j < matrixCopy.ColumnCount; j++)
                    {
                        var ww = -e[j] / e[lp1];
                        for (var ii = lp1; ii < matrixCopy.RowCount; ii++)
                        {
                            matrixCopy[ii, j] = matrixCopy[ii, j] + (ww * work[ii]);
                        }
                    }
                }

                    // Place the transformation in v for subsequent back multiplication.
                    for (i = lp1; i < matrixCopy.ColumnCount; i++)
                    {
                        vt[i, l] = e[i];
                    }
            }

            // Set up the final bidiagonal matrixCopy or order m.
            var m = Math.Min(matrixCopy.ColumnCount, matrixCopy.RowCount + 1);
            var nctp1 = nct + 1;
            var nrtp1 = nrt + 1;
            if (nct < matrixCopy.ColumnCount)
            {
                s[nctp1 - 1] = matrixCopy[(nctp1 - 1), (nctp1 - 1)];
            }

            if (matrixCopy.RowCount < m)
            {
                s[m - 1] = 0.0;
            }

            if (nrtp1 < m)
            {
                e[nrtp1 - 1] = matrixCopy[(nrtp1 - 1), (m - 1)];
            }

            e[m - 1] = 0.0;

            // If required, generate u.
           /* for (j = nctp1 - 1; j < ncu; j++)
            {
                for (i = 0; i < matrixCopy.RowCount; i++)
                {
                    u[i, j] = 0.0;
                }

                u[j, j] = 1.0;
            }

            for (l = nct - 1; l >= 0; l--)
            {
                if (s[l] != 0.0)
                {
                    for (j = l + 1; j < ncu; j++)
                    {
                        t = -DotColumns(u, matrixCopy.RowCount, l, j, l) / u[l, l];
                        for (var ii = l; ii < matrixCopy.RowCount; ii++)
                        {
                            u[ii, j] = u[ii, j] + t * u[ii, l];
                        }
                    }

                    ScalColumn(u, matrixCopy.RowCount, l, l, -1.0);
                    u[l, l]= 1.0 + u[l, l];
                    for (i = 0; i < l; i++)
                    {
                        u[i, l] = 0.0;
                    }
                }
                else
                {
                    for (i = 0; i < matrixCopy.RowCount; i++)
                    {
                        u[i, l] = 0.0;
                    }

                    u[l, l] = 1.0;
                }
            }

            // If it is required, generate v.
            for (l = matrixCopy.ColumnCount - 1; l >= 0; l--)
            {
                lp1 = l + 1;
                if (l < nrt)
                {
                    if (e[l] != 0.0)
                    {
                        for (j = lp1; j < matrixCopy.ColumnCount; j++)
                        {
                            t = -DotColumns(vt, matrixCopy.ColumnCount, l, j, lp1) / vt[lp1, l];
                            for (var ii = l; ii < matrixCopy.ColumnCount; ii++)
                            {
                                vt[ii, j] = vt[ii, j] + t * vt[ii, l];
                            }
                        }
                    }
                }

                for (i = 0; i < matrixCopy.ColumnCount; i++)
                {
                    vt[i, l] = 0.0;
                }

                vt[l, l] = 1.0;
            }*/

            // Transform s and e so that they are  double .
            /*for (i = 0; i < m; i++)
            {
                double r;
                if (s[i] != 0.0)
                {
                    t = s[i];
                    r = s[i] / t;
                    s[i] = t;
                    if (i < m - 1)
                    {
                        e[i] = e[i] / r;
                    }

                    //ScalColumn(u, matrixCopy.RowCount, i, 0, r);
                }

                // Exit
                if (i == m - 1)
                {
                    break;
                }

                if (e[i] != 0.0)
                {
                    t = e[i];
                    r = t / e[i];
                    e[i] = t;
                    s[i + 1] = s[i + 1] * r;

                    //ScalColumn(vt, matrixCopy.ColumnCount, i + 1, 0, r);
                }
            }

            // Main iteration loop for the singular values.
            /*var mn = m;
            var iter = 0;

            while (m > 0)
            {
                // Quit if all the singular values have been found. If too many iterations have been performed,
                // throw exception that Convergence Failed
                if (iter >= maxiter)
                {
                    return null;
                }

                // This section of the program inspects for negligible elements in the s and e arrays. On
                // completion the variables case and l are set as follows.
                // Case = 1     if VectorS[m] and e[l-1] are negligible and l < m
                // Case = 2     if VectorS[l] is negligible and l < m
                // Case = 3     if e[l-1] is negligible, l < m, and VectorS[l, ..., VectorS[m] are not negligible (qr step).
                // Case = 4     if e[m-1] is negligible (convergence).
                double ztest;
                double test;
                for (l = m - 2; l >= 0; l--)
                {
                    test = Math.Abs(s[l]) + Math.Abs(s[l + 1]);
                    ztest = test + Math.Abs(e[l]);
                    if (AlmostEqualRelative(ztest, test, 15))
                    {
                        e[l] = 0.0;
                        break;
                    }
                }

                int kase;
                if (l == m - 2)
                {
                    kase = 4;
                }
                else
                {
                    int ls;
                    for (ls = m - 1; ls > l; ls--)
                    {
                        test = 0.0;
                        if (ls != m - 1)
                        {
                            test = test + Math.Abs(e[ls]);
                        }

                        if (ls != l + 1)
                        {
                            test = test + Math.Abs(e[ls - 1]);
                        }

                        ztest = test + Math.Abs(s[ls]);
                        if (AlmostEqualRelative(ztest,test, 15))
                        {
                            s[ls] = 0.0;
                            break;
                        }
                    }

                    if (ls == l)
                    {
                        kase = 3;
                    }
                    else if (ls == m - 1)
                    {
                        kase = 1;
                    }
                    else
                    {
                        kase = 2;
                        l = ls;
                    }
                }

                l = l + 1;

                // Perform the task indicated by case.
                /*int k;
                double f;
                double sn;
                double cs;
                switch (kase)
                {
                    // Deflate negligible VectorS[m].
                    case 1:
                        f = e[m - 2];
                        e[m - 2] = 0.0;
                        double t1;
                        for (var kk = l; kk < m - 1; kk++)
                        {
                            k = m - 2 - kk + l;
                            t1 = s[k];
                            Rotg(ref t1, ref f, out cs, out sn);
                            s[k] = t1;
                            if (k != l)
                            {
                                f = -sn * e[k - 1];
                                e[k - 1] = cs * e[k - 1];
                            }

                            //RotColumns(vt, matrixCopy.ColumnCount, k, m - 1, cs, sn);
                        }

                        break;

                    // Split at negligible VectorS[l].
                    /*case 2:
                        f = e[l - 1];
                        e[l - 1] = 0.0;
                        for (k = l; k < m; k++)
                        {
                            t1 = s[k];
                            Rotg(ref t1, ref f, out cs, out sn);
                            s[k] = t1;
                            f = -sn * e[k];
                            e[k] = cs * e[k];

                           // RotColumns(u, matrixCopy.RowCount, k, l - 1, cs, sn);
                        }

                        break;

                    // Perform one qr step.
                    /*case 3:
                        // Calculate the shift.
                        var scale = 0.0;
                        scale = Math.Max(scale, Math.Abs(s[m - 1]));
                        scale = Math.Max(scale, Math.Abs(s[m - 2]));
                        scale = Math.Max(scale, Math.Abs(e[m - 2]));
                        scale = Math.Max(scale, Math.Abs(s[l]));
                        scale = Math.Max(scale, Math.Abs(e[l]));
                        var sm = s[m - 1] / scale;
                        var smm1 = s[m - 2] / scale;
                        var emm1 = e[m - 2] / scale;
                        var sl = s[l] / scale;
                        var el = e[l] / scale;
                        var b = (((smm1 + sm) * (smm1 - sm)) + (emm1 * emm1)) / 2.0;
                        var c = (sm * emm1) * (sm * emm1);
                        var shift = 0.0;
                        if (b != 0.0 || c != 0.0)
                        {
                            shift = Math.Sqrt((b * b) + c);
                            if (b < 0.0)
                            {
                                shift = -shift;
                            }

                            shift = c / (b + shift);
                        }

                        f = ((sl + sm) * (sl - sm)) + shift;
                        var g = sl * el;

                        // Chase zeros.
                        /*for (k = l; k < m - 1; k++)
                        {
                            Rotg(ref f, ref g, out cs, out sn);
                            if (k != l)
                            {
                                e[k - 1] = f;
                            }

                            f = (cs * s[k]) + (sn * e[k]);
                            e[k] = (cs * e[k]) - (sn * s[k]);
                            g = sn * s[k + 1];
                            s[k + 1] = cs * s[k + 1];

                            //RotColumns(vt, matrixCopy.ColumnCount, k, k + 1, cs, sn);

                           /* Rotg(ref f, ref g, out cs, out sn);
                            s[k] = f;
                            f = (cs * e[k]) + (sn * s[k + 1]);
                            s[k + 1] = (-sn * e[k]) + (cs * s[k + 1]);
                            g = sn * e[k + 1];
                            e[k + 1] = cs * e[k + 1];
                            if (k < matrixCopy.RowCount)
                            {
                                RotColumns(u, matrixCopy.RowCount, k, k + 1, cs, sn);
                            }
                        }

                        e[m - 2] = f;
                        iter = iter + 1;
                        break;

                    // Convergence.
                   /* case 4:
                        // Make the singular value  positive
                        /*if (s[l] < 0.0)
                        {
                            s[l] = -s[l];

                            //ScalColumn(vt, matrixCopy.ColumnCount, l, 0, -1.0);
                        }

                        // Order the singular value.
                       /* while (l != mn - 1)
                        {
                            if (s[l] >= s[l + 1])
                            {
                                break;
                            }

                            t = s[l];
                            s[l] = s[l + 1];
                            s[l + 1] = t;
                            if (l < matrixCopy.ColumnCount)
                            {
                                SwapColumns(vt, matrixCopy.ColumnCount, l, l + 1);
                            }

                            if (l < matrixCopy.RowCount)
                            {
                                SwapColumns(u, matrixCopy.RowCount, l, l + 1);
                            }

                            l = l + 1;
                        }

                        iter = 0;
                        m = m - 1;
                        break;
                }
            }

            //vt = vt.Transpose();

            // Adjust the size of s if rows < columns. We are using ported copy of linpack's svd code and it uses
            // a singular vector of length mRows+1 when mRows < mColumns. The last element is not used and needs to be removed.
            // we should port lapack's svd routine to remove this problem.
            /*if (matrixCopy.RowCount < matrixCopy.ColumnCount)
            {
                nm--;
                var tmp = new Vector(nm);
                for (i = 0; i < nm; i++)
                {
                    tmp[i] = s[i];
                }

                s = tmp;
            }

            return new SVD(s, u, vt);
        }

        static double Dsign(double z1, double z2)
        {
            return Math.Abs(z1) * (z2 / Math.Abs(z2));
        }

        static void SwapColumns(Matrix a, int rowCount, int columnA, int columnB)
        {
            for (var i = 0; i < rowCount; i++)
            {
                double z = a[i, columnA];
                a[i, columnA] = a[i, columnB];
                a[i, columnB] = z;
            }
        }

        static double NrmColumn(Matrix a, int rowCount, int column, int rowStart)
        {
            double s = 0;
            for (var i = rowStart; i < rowCount; i++)
            {
                s += a[i, column] * a[i, column];
            }

            return Math.Sqrt(s);
        }

        static double DotColumns(Matrix a, int rowCount, int columnA, int columnB, int rowStart)
        {
            var z = 0.0;
            for (var i = rowStart; i < rowCount; i++)
            {
                z += a[i, columnB] * a[i, columnA];
            }

            return z;
        }

        static void RotColumns(Matrix a, int rowCount, int columnA, int columnB, double c, double s)
        {
            for (var i = 0; i < rowCount; i++)
            {
                var z = (c * a[i, columnA]) + (s * a[i, columnB]);
                var tmp = (c * a[i, columnB]) - (s * a[i, columnA]);
                a[i, columnB] = tmp;
                a[i, columnA] = z;
            }
        }

        static void ScalColumn(Matrix a, int rowCount, int column, int rowStart, double z)
        {
            for (var i = rowStart; i < rowCount; i++)
            {
                a[i, column] = a[i, column] * z;
            }
        }

        /*public static bool AlmostEqualRelative(double a, double b, double maximumError)
        {
            return AlmostEqualNormRelative(a, b, a - b, maximumError);
        }

        /*static bool AlmostEqualNormRelative(double a, double b, double diff, double maximumError)
        {
            // If A or B are infinity (positive or negative) then
            // only return true if they are exactly equal to each other -
            // that is, if they are both infinities of the same sign.
            /*if (double.IsInfinity(a) || double.IsInfinity(b))
            {
             *   return a == b;
            }*/

            // If A or B are a NAN, return false. NANs are equal to nothing,
            // not even themselves.
            /*if (double.IsNaN(a) || double.IsNaN(b))
            {
              * return false;
            }*/

            // If one is almost zero, fall back to absolute equality
            /*if (Math.Abs(a) < Math.Pow(2, -53) || Math.Abs(b) < Math.Pow(2, -53))
            {
                return Math.Abs(diff) < maximumError;
            }

            if ((a == 0 && Math.Abs(b) < maximumError) || (b == 0 && Math.Abs(a) < maximumError))
            {
                return true;
            }

            return Math.Abs(diff) < maximumError * Math.Max(Math.Abs(a), Math.Abs(b));
        }*/
        /*static void Rotg(ref double da, ref double db, out double c, out double s)
        {
            double r, z;

            var roe = db;
            var absda = Math.Abs(da);
            var absdb = Math.Abs(db);
            if (absda > absdb)
            {
                roe = da;
            }

            var scale = absda + absdb;
            if (scale == 0.0)
            {
                c = 1.0;
                s = 0.0;
                r = 0.0;
                z = 0.0;
            }
            else
            {
                var sda = da / scale;
                var sdb = db / scale;
                r = scale * Math.Sqrt((sda * sda) + (sdb * sdb));
                if (roe < 0.0)
                {
                    r = -r;
                }

                c = da / r;
                s = db / r;
                z = 1.0;
                if (absda > absdb)
                {
                    z = s;
                }

                if (absdb >= absda && c != 0.0)
                {
                    z = 1.0 / c;
                }
            }

            da = r;
            db = z;
        }*/
    }
}
