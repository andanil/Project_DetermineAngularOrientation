using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Maths.LinearAlgebra
{
    public class MatrixOperations
    {
        public static Eigendecomp Eigendecomposition(Matrix matrix, double eps)
        {
            Eigendecomp result = new Eigendecomp();
            result.Compute(matrix, eps);
            return result;
        }
    }
}
