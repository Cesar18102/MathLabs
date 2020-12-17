using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;

using ScottPlot;

namespace MathLab4
{
    [Obsolete]
    public class Program
    {
        public static void Main(string[] args)
        {
            Program program = new Program();

            Form form = new Form();

            form.AutoSize = true;
            form.AutoScroll = true;

            program.DoTasks(form);

            form.ShowDialog();
        }

        public void DoTasks(Form form)
        {
            Point point1 = Task1(form, new Point(0, 0));
        }

        public Point Task1(Form form, Point location)
        {
            double a = 0.9;
            double[] f = new double[] { 4 * a, -4.4 * a, -1, 1 };

            double[] firstDerivative = Poly.GetPolyDerivative(f);
            double[] secondDerivative = Poly.GetPolyDerivative(firstDerivative);

            double[] xs = Range.Make(-5, 5, 0.1);
            double[] breakXs = Poly.PolySolvePrimitive(secondDerivative, 0);

            (double[] mins, double[] maxs) = FindExtremumsUniform(f, -4, 4, 0.01);

            Console.WriteLine("MINS: " + string.Join(" ", mins));
            Console.WriteLine("MAXS: " + string.Join(" ", maxs));

            Plot plot = new Plot();

            plot.PlotSignalXY(xs, Poly.Val(f, xs), Color.Red, 7, label: "Original signal");
            plot.PlotSignalXY(xs, Poly.Val(firstDerivative, xs), Color.Green, 5, label: "First Derivative");
            plot.PlotSignalXY(xs, Poly.Val(secondDerivative, xs), Color.Blue, 3, label: "Second Derivative");

            plot.PlotScatter(
                breakXs.Concat(breakXs).Concat(breakXs).ToArray(),
                Poly.Val(f, breakXs).Concat(Poly.Val(firstDerivative, breakXs)).Concat(Poly.Val(secondDerivative, breakXs)).ToArray(),
                Color.Yellow, label: "Unimodal break points", markerSize: 3, lineStyle: LineStyle.Dot
            );

            plot.PlotScatter(mins, Poly.Val(f, mins), Color.Cyan, markerSize: 7, label: "Min points", lineStyle: LineStyle.Dot);
            plot.PlotScatter(maxs, Poly.Val(f, maxs), Color.LightGreen, markerSize: 7, label: "Max points", lineStyle: LineStyle.Dot);

            plot.Legend();

            return plot.AddToForm(form, location);
        }

        public (double[] mins, double[] maxs) FindExtremumsUniform(double[] poly, double start, double end, double delta)
        {
            List<double> mins = new List<double>();
            List<double> maxs = new List<double>();

            double lastX = start;
            while(lastX <= end)
            {
                Extremum extremum = FindExtremumUniformForUnimodal(poly, lastX, end, delta);

                if (extremum == null)
                    break;

                if (extremum.IsMin)
                    mins.Add(extremum.X);
                else
                    maxs.Add(extremum.X);

                lastX = extremum.X;
            }

            return (mins.ToArray(), maxs.ToArray());
        }

        public Extremum FindExtremumUniformForUnimodal(double[] poly, double start, double end, double delta)
        {
            double[] xs = Range.Make(start, end, delta);
            Extremum extremum = new Extremum() { IsMin = true, X = start };

            Predicate<int> breakCondition = i => Poly.Val(poly, xs[i + 1])[0] >= Poly.Val(poly, xs[i])[0];
            if (breakCondition(0))
            {
                breakCondition = i => Poly.Val(poly, xs[i + 1])[0] <= Poly.Val(poly, xs[i])[0];
                extremum.IsMin = false;
            }

            for (int i = 0; i < xs.Length - 1; ++i)
                if (breakCondition(i))
                {
                    extremum.X = xs[i] + delta / 2;
                    return extremum;
                }

            return null;
        }

        public class Extremum
        {
            public double X { get; set; }
            public bool IsMin { get; set; }
        }

        public static class Poly
        {
            public static double[] PolySolvePrimitive(double[] poly, double eq)
            {
                if (poly.Length == 2)
                    return new double[] { (eq - poly[0]) / poly[1] };

                throw new ArgumentException();
            }

            public static double[] GetPolyDerivative(double[] poly)
            {
                double[] diffPoly = new double[poly.Length - 1];

                for (int i = 0; i < diffPoly.Length; ++i)
                    diffPoly[i] = poly[i + 1] * (i + 1);

                return diffPoly;
            }

            public static Func<double, double> Get(double[] p)
            {
                return x => p.Select((a, i) => a * Math.Pow(x, i)).Sum();
            }

            public static double[] Val(double[] p, params double[] x)
            {
                Func<double, double> poly = Get(p);
                return x.Select(xi => poly(xi)).ToArray();
            }
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

    public static class Extensions
    {
        public static Point AddToForm(this Plot plot, Form form, Point location)
        {
            Bitmap chart = plot.Render();

            PictureBox picture = new PictureBox();

            picture.Size = chart.Size;
            picture.Location = location;
            picture.BackgroundImage = chart;

            form.Controls.Add(picture);

            return new Point(location.X, location.Y + chart.Height);
        }
    }
}
