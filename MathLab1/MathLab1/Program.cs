using System;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;

using ScottPlot;

namespace MathLab1
{
    public class Program
    {
        public static void Main(string[] args)
        {
            Program program = new Program();

            program.Task1();
            program.Task2();
            program.Task3();
            program.Task4();
        }

        public (int height, int width) GetMatrixSize<T>(T[,] matrix)
        {
            return (matrix.GetLength(0), matrix.GetLength(1));
        }

        public void PrintMatrix<T>(T[,] matrix)
        {
            var size = GetMatrixSize(matrix);
            for(int i = 0; i < size.height; ++i)
            {
                for (int j = 0; j < size.width; ++j)
                    Console.Write($"{matrix[i, j]} ");
                Console.WriteLine();
            }
        }

        public void PrintVector<T>(T[] vector, string name)
        {
            Console.WriteLine($"{name} = {string.Join(" ", vector)}");
        }

        public double[,] ConvertToDouble(int[,] matrix)
        {
            var size = GetMatrixSize(matrix);
            double[,] result = new double[size.height, size.width];
            for (int i = 0; i < size.height; ++i)
                for (int j = 0; j < size.width; ++j)
                    result[i, j] = (double)matrix[i, j];
            return result;
        }

        public double[] ConvertToDouble(int[] vector)
        {
            return vector.Select(item => (double)item).ToArray();
        }

        public T[,] VectorToMatrix<T>(T[] vector, bool inline)
        {
            if(inline)
            {
                T[,] result = new T[1, vector.Length];
                for (int i = 0; i < vector.Length; ++i)
                    result[0, i] = vector[i];
                return result;
            }
            else
            {
                T[,] result = new T[vector.Length, 1];
                for (int i = 0; i < vector.Length; ++i)
                    result[i, 0] = vector[i];
                return result;
            }
        }

        public double[,] Multiply(double[,] a, double[,] b)
        {
            var sizea = GetMatrixSize(a);
            var sizeb = GetMatrixSize(b);

            if (sizea.width != sizeb.height)
                throw new ArgumentException();

            double[,] result = new double[sizea.height, sizeb.width];

            for(int i = 0; i < sizea.height; ++i)
                for(int j = 0; j < sizeb.width; ++j)
                    for(int k = 0; k < sizea.width; ++k)
                        result[i, j] += a[i, k] * b[k, j];

            return result;
        }

        public void Task1()
        {
            Console.WriteLine($"7 + 6 = {7 + 6}");
            Console.WriteLine($"7 * 6 = {7 * 6}");

            Func<double, double> f = new Func<double, double>(x => x * x + 3 * x - 1);
            Console.WriteLine($"\nf(x) = x^2 + 3x - 1;\nf(5) = {f(5)}");
        }

        public void Task2()
        {
            int count = 4;

            int[] k = Enumerable.Range(1, count).ToArray();
            PrintVector(k, "\nK");

            int[] p = new int[] { 3, 5, 1, 4 };
            PrintVector(p, "P");

            int[] sum = Enumerable.Range(0, count).Select(i => k[i] + p[i]).ToArray();
            PrintVector(sum, "K + P");

            int mult = Enumerable.Range(0, count).Aggregate(
                0, (a, i) => a + k[i] * p[i]
            );
            Console.WriteLine($"K * P = {mult}\n");

            int[,] m = new int[4, 4] 
            {
                { 1, 3, 2, 7 },
                { 5, 6, 3, 2 },
                { 8, 9, 11, 1 },
                { -1, 7, 3, 0 }
            };
            Console.WriteLine("M = ");
            PrintMatrix(m);

            double[,] result = Multiply(
                ConvertToDouble(m),
                ConvertToDouble(VectorToMatrix(p, false))
            );
            Console.WriteLine("\nM * P = ");
            PrintMatrix(result);
        }

        public T[] MatrixToVector<T>(T[,] matrix)
        {
            var size = GetMatrixSize(matrix);

            if(size.height == 1)
            {
                T[] vector = new T[size.width];
                for (int i = 0; i < size.width; ++i)
                    vector[i] = matrix[0, i];
                return vector;
            }
            else if(size.width == 1)
            {
                T[] vector = new T[size.height];
                for (int i = 0; i < size.height; ++i)
                    vector[i] = matrix[i, 0];
                return vector;
            }

            throw new ArgumentException();
        }

        public T[,] GetMinor<T>(T[,] matrix, int row, int column)
        {
            var size = GetMatrixSize(matrix);
            T[,] submatrix = new T[size.height - 1, size.width - 1];

            int di = 0;
            for(int i = 0; i < size.height; ++i)
            {
                if(i == row)
                {
                    di = -1;
                    continue;
                }

                int dj = 0;
                for(int j = 0; j < size.width; ++j)
                {
                    if(j == column)
                    {
                        dj = -1;
                        continue;
                    }

                    submatrix[i + di, j + dj] = matrix[i, j];
                }
            }

            return submatrix;
        }

        public double Det(double[,] matrix)
        {
            var size = GetMatrixSize(matrix);

            if (size.height != size.width)
                throw new ArgumentException();

            if (size.width == 2 && size.height == 2)
                return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];

            return Enumerable.Range(0, size.width).Sum(
                i => Math.Pow(-1, i) * matrix[0, i] * Det(GetMinor(matrix, 0, i))
            );
        }

        public T[,] Transpose<T>(T[,] matrix)
        {
            var size = GetMatrixSize(matrix);
            T[,] transposed = new T[size.width, size.height];
            for (int i = 0; i < size.height; ++i)
                for (int j = 0; j < size.width; ++j)
                    transposed[j, i] = matrix[i, j];
            return transposed;
        }

        public void Multiply(double[,] matrix, double value)
        {
            var size = GetMatrixSize(matrix);
            for (int i = 0; i < size.height; ++i)
                for (int j = 0; j < size.width; ++j)
                    matrix[i, j] *= value;
        }

        public double[,] Invert(double[,] matrix)
        {
            double det = Det(matrix);

            if (det == 0)
                throw new ArgumentException();

            var size = GetMatrixSize(matrix);
            double[,] transposed = Transpose(matrix);
            double[,] complements = new double[size.height, size.width];

            for (int i = 0; i < size.height; ++i)
                for (int j = 0; j < size.width; ++j)
                    complements[i, j] = Math.Pow(-1, i + j) * Det(GetMinor(transposed, i, j));

            Multiply(complements, 1 / det);

            return complements;
        }

        public double[] PolyFit(double[] x, double[] y, int k)
        {
            if (x.Length != y.Length)
                throw new ArgumentException();

            int n = x.Length;

            double[] powersums = new double[2 * k + 1];
            for (int i = 0; i <= 2 * k; ++i)
                powersums[i] = Enumerable.Range(0, n).Sum(j => Math.Pow(x[j], i));

            double[,] xm = new double[k + 1, k + 1];
            double[] ym = new double[k + 1];

            for (int i = 0; i <= k; ++i)
            {
                for (int j = 0; j <= k; ++j)
                    xm[i, j] = powersums[i + j];
                ym[i] = Enumerable.Range(0, n).Sum(j => Math.Pow(x[j], i) * y[j]);
            }

            return MatrixToVector(
                Multiply(Invert(xm), VectorToMatrix(ym, false))
            );
        }

        public Func<double, double> Poly(double[] p)
        {
            return x => p.Select((a, i) => a * Math.Pow(x, i)).Sum();
        }

        public double[] PolyVal(double[] p, double[] x)
        {
            Func<double, double> poly = Poly(p);
            return x.Select(xi => poly(xi)).ToArray();
        }

        public void OutPolyFitVal(double[] x, double[] y, int n, Plot plot, Plot plotSmooth, double smooth, Color color, double thickness = 1)
        {
            double[] pfn = PolyFit(x, y, n);
            double[] pvn = PolyVal(pfn, x);

            PrintVector(pfn, $"\nQ{n}");
            PrintVector(pvn, $"V{n}");

            plot.PlotSignalXY(x, pvn, color, thickness);

            double min = x.Min();
            double spread = x.Max() - min;
            int count = (int)(spread / smooth);
            double step = spread / count;

            double[] xsmooth = Enumerable.Range(0, count).Select(i => i * step + min).ToArray();
            double[] ysmooth = PolyVal(pfn, xsmooth);

            plotSmooth.PlotSignalXY(xsmooth, ysmooth, color, thickness);
        }

        public void ShowPlotForm(params Bitmap[] renderedPlots)
        {
            Form plotForm = new Form();

            plotForm.AutoScroll = true;
            plotForm.Size = new Size(
                renderedPlots.Take(2).Sum(plot => plot.Width),
                renderedPlots.Max(plot => plot.Height)
            );

            Point location = new Point(0, 0);
            foreach(Bitmap plot in renderedPlots)
            {
                PictureBox picture = new PictureBox();

                picture.Size = plot.Size;
                picture.Location = location;
                picture.Image = plot;

                plotForm.Controls.Add(picture);
                location.Offset(plot.Width, 0);
            }

            plotForm.ShowDialog();
        }

        public void Task3()
        {
            Func<double, double> y1 = xv => Math.Sin(xv);
            Func<double, double> y2 = xv => 1.0 / xv;

            double[] xs = Range.Make(0.1, 1, 0.05);

            double[] y1s = xs.Select(xv => y1(xv)).ToArray();
            double[] y2s = xs.Select(xv => y2(xv)).ToArray();

            Plot plotY1 = new Plot();
            plotY1.PlotSignalXY(xs, y1s);

            Plot plotY2 = new Plot();
            plotY2.PlotSignalXY(xs, y2s);

            ShowPlotForm(
                plotY1.Render(),
                plotY2.Render()
            );
        }

        public void Task4()
        { 
            Plot plotPoints = new Plot();
            Plot plotSmooth = new Plot();

            double[] x = new double[] { -2.1, -1, 0.5, 1, 2.1, 3.1 };
            double[] y = new double[] { 0.6, 0.1, -0.5, 0.1, 0.7, 0.9 };

            PrintVector(x, $"\nX");
            PrintVector(y, $"Y");

            plotPoints.PlotSignalXY(x, y, Color.Violet, 5);
            plotSmooth.PlotSignalXY(x, y, Color.Violet, 5);

            OutPolyFitVal(x, y, 2, plotPoints, plotSmooth, 0.05, Color.Green);
            OutPolyFitVal(x, y, 3, plotPoints, plotSmooth, 0.05, Color.Red);
            OutPolyFitVal(x, y, 5, plotPoints, plotSmooth, 0.05, Color.Blue);

            ShowPlotForm(
                plotPoints.Render(), 
                plotSmooth.Render()
            );
        }

        public static class Range
        {
            public static double[] Make(double start, double end, double delta)
            {
                int count = (int)Math.Ceiling((end - start) / delta);
                return Enumerable.Range(0, count).Select(i => i * delta + start).ToArray();
            }
        }
    }
}
