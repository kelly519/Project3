using System;

namespace InterpolationPolynomial
{
    class Program
    {
        static void Main(string[] args)
        {
            int i = 0;
            List<Tuple<double, double>> points = new List<Tuple<double, double>>();

            Console.WriteLine("Enter the points in the format x, y. Enter END to stop.");

            while (true)
            {
                ++i;
                Console.Write("Point #{0}: ", i);
                string input = Console.ReadLine();
                if (input.ToUpper() == "END")
                {
                    break;
                }

                string[] split = input.Split(' ');
                double x = double.Parse(split[0]);
                double y = double.Parse(split[1]);
                points.Add(new Tuple<double, double>(x, y));
            }

            int n = points.Count - 1;
            Console.WriteLine("The order of the interpolation polynomial is " + n);

            double[,] A = new double[n + 1, n + 1];
            double[] B = new double[n + 1];

            for (int k = 0; k <= n; k++)
            {
                B[k] = points[k].Item2;
                for (int j = 0; j <= n; j++)
                {
                    A[k, j] = Math.Pow(points[k].Item1, j);
                }
            }

            double[] coefficients = GaussianElimination(A, B);

            Console.WriteLine("The obtained polynomial is:");
            for (int k = n; k >= 0; k--)
            {
                Console.Write(coefficients[k].ToString("0.00000") + "x^" + k + " ");
                if (k > 0)
                {
                    Console.Write("+ ");
                }
            }

            Console.WriteLine("\nValues of the polynomial for -1, 0 and 1:");
            Console.WriteLine("f(1) = " + EvaluatePolynomial(coefficients, 1.0));
            Console.WriteLine("f(0) = " + EvaluatePolynomial(coefficients, 0.0));
            Console.WriteLine("f(-1) = " + EvaluatePolynomial(coefficients, -1.0));

            Console.WriteLine("\nThe derivative of the polynomial is:");

            double[] derivativeCoefficients = CalculateDerivative(coefficients);
            for (int k = derivativeCoefficients.Length - 1; k >= 1; k--)
            {
                Console.Write(derivativeCoefficients[k].ToString("0.00000") + "x^" + (k - 1) + " + ");
            }
            Console.WriteLine(derivativeCoefficients[0].ToString("0.00000"));
            Console.WriteLine("\nLooking for a root with inital guess 2");
            double guess = 2;
            double tolerance = 1e-6;
            int maxIterations = 100;

            double root = NewtonMethod(coefficients, derivativeCoefficients, guess, tolerance, maxIterations);
            Console.WriteLine("A root of the polynomial is: " + root.ToString("0.00000"));
        }

        static double NewtonMethod(double[] coefficients, double[] derivativeCoefficients, double guess, double tolerance, int maxIterations)
        {
            double x = guess;
            int iteration = 0;
            double error = double.MaxValue;
            while (error > tolerance && iteration < maxIterations)
            {
                double fx = EvaluatePolynomial(coefficients, x);
                double fDerivative = EvaluatePolynomial(derivativeCoefficients, x);
                double xNew = x - fx / fDerivative;
                error = Math.Abs(xNew - x);
                x = xNew;
                iteration++;
            }
            return x;
        }

        static double[] CalculateDerivative(double[] coefficients)
        {
            int n = coefficients.Length - 1;
            double[] derivativeCoefficients = new double[n];
            for (int i = 1; i <= n; i++)
            {
                derivativeCoefficients[i - 1] = coefficients[i] * i;
            }
            return derivativeCoefficients;
        }

        static double[] GaussianElimination(double[,] A, double[] B)
        {
            int n = B.Length;
            for (int i = 0; i < n; i++)
            {
                int max = i;
                for (int j = i + 1; j < n; j++)
                {
                    if (Math.Abs(A[j, i]) > Math.Abs(A[max, i]))
                    {
                        max = j;
                    }
                }

                double[] temp = new double[n + 1];
                for (int j = 0; j < n; j++)
                {
                    temp[j] = A[i, j];
                    A[i, j] = A[max, j];
                    A[max, j] = temp[j];
                }

                double t = B[i];
                B[i] = B[max];
                B[max] = t;

                for (int j = i + 1; j < n; j++)
                {
                    double factor = A[j, i] / A[i, i];
                    B[j] -= factor * B[i];
                    for (int k = i; k < n; k++)
                    {
                        A[j, k] -= factor * A[i, k];
                    }
                }
            }

            double[] x = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = B[i];
                for (int j = i + 1; j < n; j++)
                {
                    x[i] -= A[i, j] * x[j];
                }
                x[i] = x[i] / A[i, i];
            }

            return x;
        }

        static double EvaluatePolynomial(double[] coefficients, double x)
        {
            double result = 0;
            for (int i = 0; i < coefficients.Length; i++)
            {
                result += coefficients[i] * Math.Pow(x, i);
            }
            return result;
        }
    }
}