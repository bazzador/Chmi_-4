using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using MathNet.Numerics.LinearAlgebra.Double;
using OxyPlot;
using OxyPlot.Series;
using System.Collections;
using System.Diagnostics;
using System.Threading;
using Antlr.Runtime;
using NCalc.Domain; 
using OxyPlot.WindowsForms;
using System.Windows.Controls.Primitives;
namespace WindowsFormsApp1
{
    public partial class Form1 : Form
    {
        private static double[] x = { 2, 5, 7, 9, 12 };
        private static double[] y = { 2, 0, -1, 1, 10 };
        public double[] h = new double[4];
        public double[] u;
        public double[] d = new double[4];
        public List<string> expressions;
        private double[] spl, lagr;

        public Form1()
        {
            InitializeComponent();
            InitializePlot();
        }
        private void InitializePlot()
        {
            var model = new PlotModel { Title = "Сплайн-функція" };
            plotView2 = new PlotView
            {
                Dock = DockStyle.Fill
            };
            Controls.Add(plotView2);
            foreach (var series in SLAR())
            {
                model.Series.Add(series);
            }
            plotView2.Model = model;
            plotView2 = new PlotView
            {
                Dock = DockStyle.Fill
            };
            Controls.Add(plotView2);
            var lagrangeSeries = GenerateLagrangeSeries();
            model.Series.Add(lagrangeSeries);
        }
            public void hud()
            {

                for (int k = 0; k < x.Length - 1; k++)
                {
                    h[k] = x[k + 1] - x[k];
                    d[k] = (y[k + 1] - y[k]) / h[k];
                }
                 u = new double[d.Length - 1];
                for (int k = 1; k < u.Length + 1; k++)
                {
                    u[k - 1] = 6 * (d[k] - d[k - 1]);
                }
            }

        private IEnumerable<FunctionSeries> SLAR()
            {
                hud();
                double[,] slar = new double[,]
                {
                    { 2*(h[0]+h[1]), h[1], 0},
                    { h[1], 2*(h[1]+h[2]), h[2]},
                    { 0, h[2], 2*(h[2]+h[3])}
                };

                double[] b = new double[] {u[0], u[1], u[2]};
                double[] m1 = Solve(slar, b);
                double[] m = new double[5];
                m[1] = m1[0];
                m[2] = m1[1];
                m[3] = m1[2];
                string Sx = $"Кубічний сплайн\n";
                expressions = new List<string>();
                for (int i = 0; i < 4; i++)
                {
                      Sx = ($"{(m[i+1] - m[i]) / (6 * h[i])}*Pow((x - {x[i]}),3)  + {m[i]/2}*Pow((x - {x[i]}),2) + {(d[i] - (h[i] * (2 * m[i] + m[i+1]))/6)}*(x - {x[i]}) + {y[i]}\n").ToString();
                      Sx = Sx.Replace(',', '.');
                      expressions.Add(Sx);
                }
                int c = 1;
                int q = 0;
                spl = new double[0];
                foreach (var expression in expressions)
                {
                    var series = new FunctionSeries();
                    var exp = new NCalc.Expression(expression.Replace(").", "),"));
                    double startX = 0;
                    double endX = 0;
                    switch (c)
                    {
                        case 1:
                            startX = x[0];
                            endX = x[1];
                            break;
                        case 2:
                            startX = x[1];
                            endX = x[2];
                            break;
                        case 3:
                            startX = x[2];
                            endX = x[3];
                            break;
                        case 4:
                            startX = x[3];
                            endX = x[4];
                            break;
                        default:
                            break;
                    }
                    double y = 0;
                    for (double x = startX; ; x += 0.1)
                    {
                        x = Math.Round(x, 1);
                        if (x > endX)
                        {
                            break;
                        }
                        Array.Resize(ref spl, q + 1);
                        exp.Parameters["x"] = x;
                        y = Convert.ToDouble(exp.Evaluate());
                        series.Points.Add(new DataPoint(x, y));
                        spl[q] = y;
                        q++;
                    }
                    c++;
                    yield return series;
                }
            }
        private FunctionSeries GenerateLagrangeSeries()
        {
            lagr = new double[0];
            double[] Y = { 1, 14, 8, 12, 10 };
            int w = 0;
            var lagrangePoints = new List<DataPoint>();
            for (double x_ = x[0]; ; x_ += 0.1)
            {
                x_ = Math.Round(x_, 1);
                if (x_ > x[x.Length - 1])
                {
                    break;
                }
                Array.Resize(ref lagr, w + 1);
                double y = LagrangeInterpolation(x_, x, Y);

                lagrangePoints.Add(new DataPoint(x_, y));
                lagr[w] = y;
                w++;

            }
            var lagrangeSeries = new FunctionSeries();
            foreach (var point in lagrangePoints)
            {
                lagrangeSeries.Points.Add(point);
            }
            return lagrangeSeries;
        }

    private double LagrangeInterpolation(double x, double[] xValues, double[] yValues)
    {
        double result = 0;
        for (int i = 0; i < xValues.Length; i++)
        {
            double term = yValues[i];
            for (int j = 0; j < xValues.Length; j++)
            {
                if (j != i)
                    term *= (x - xValues[j]) / (xValues[i] - xValues[j]);
            }
            result += term;
        }
        return result;
    }
        public double[] Solve(double[,]slar, double[]b)
            {
                var matrixA = DenseMatrix.OfArray(slar);
                var vectorB = DenseVector.OfArray(b);

                var lu = matrixA.LU();
                var solution = lu.Solve(vectorB);

                MathNet.Numerics.LinearAlgebra.Storage.VectorStorage<double> vectorStorage = solution.Storage;

                int line = vectorStorage.Length;
                double[] result = new double[line];

                for (int i = 0; i < result.Length;i++)
                {
                    result[i] = vectorStorage[i];
                }
                return result;
            }
        private FunctionSeries Diff()
        {
            int c = 0;
            var diff = new List<DataPoint>();
            for (double i = x[0]; ; i += 0.2)
            {
                i = Math.Round(i, 1);
                if (i > x[x.Length - 1])
                {
                    break;
                }
                diff.Add(new DataPoint(i, Math.Abs(lagr[c] - spl[c])));
                c++;
            }
            var diffSeries = new FunctionSeries();
            foreach (var point in diff)
            {
                diffSeries.Points.Add(point);
            }
            return diffSeries;
        }
        private void label1_Click(object sender, EventArgs e)
        {
            Form2 form2 = new Form2();
            var model = new PlotModel { Title = "Різниця" };
            PlotView plotView = new PlotView
            {
                Dock = DockStyle.Fill
            };
            form2.Controls.Add(plotView);
            var dif = Diff();
            model.Series.Add(dif);
            plotView.Model = model;
            form2.textBox1.Visible = false;
            form2.ShowDialog();
        }
        private void label2_Click(object sender, EventArgs e)
        {
            Form2 form2 = new Form2();
            form2.textBox1.AppendText($"Lagrange formula" + Environment.NewLine);
            form2.textBox1.AppendText($"Spline formula" + Environment.NewLine);
            form2.textBox1.AppendText($"{expressions[0].Replace(").", "),")}" + Environment.NewLine);
            form2.textBox1.AppendText(Environment.NewLine);
            form2.textBox1.AppendText($"{expressions[1].Replace(").", "),")}" + Environment.NewLine);
            form2.textBox1.AppendText(Environment.NewLine);
            form2.textBox1.AppendText($"{expressions[2].Replace(").", "),")}" + Environment.NewLine);
            form2.textBox1.AppendText(Environment.NewLine);
            form2.textBox1.AppendText($"{expressions[3].Replace(").", "),")}" + Environment.NewLine);
            form2.textBox1.ReadOnly = true;
            form2.ShowDialog();
        }
    }
}
