using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Numerics;
using System.Reflection;
using System.Reflection.Metadata;
using System.Runtime.ExceptionServices;
using System.Runtime.InteropServices.Marshalling;
using System.Runtime.Intrinsics.X86;
using System.Text;
using static System.Runtime.InteropServices.JavaScript.JSType;

public class Program
{
    public static void Main()
    {
        Program program = new Program();
    }

    void FindMaxIndex(int[,] matrix, out int row, out int column)
    {
        row = 0; column = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
                if (matrix[i, j] > matrix[row, column]) { row = i; column = j; }
    }
    int FindMaxInColumn(int[,] matr)
    {
        int max = -1000000, stroka = 0;
        for (int i = 0; i < matr.GetLength(0); i++)
        {
            if (matr[i, 0] > max)
            {
                max = matr[i, 0];
                stroka = i;
            }
        }
        return stroka;
    }
    static int FindMaxInStolb(int[,] matr)
    {
        int max = -1000000, stolb = 0;
        for (int j = 0; j < matr.GetLength(1); j++)
        {
            if (matr[0, j] > max)
            {
                max = matr[0, j];
                stolb = j;
            }
        }
        return stolb;
    }
    #region Level 1
    public long Task_1_1(int n, int k)
    {
        long answer = 0;
        static long Factorial(int n)
        {
            long output = 1;
            for (int i = 1; i < n + 1; i++) output *= i;
            return output;
        }
        static long Combinations(int n, int k)
        {
            long output = 1;
            output *= Factorial(n) / (Factorial(k) * Factorial(n - k));
            return output;
        }
        if (k < 0 || k > n)
        {
            return answer = 0;
        }
        answer = Combinations(n, k);
        return answer;
    }

    public int Task_1_2(double[] first, double[] second)
    {
        int answer = 0;
        static double GeronArea(double a, double b, double c)
        {
            double p = (a + b + c) / 2;
            return Math.Sqrt(p * (p - a) * (p - b) * (p - c));
        }
        if ((first[0] >= first[1] + first[2]) || (first[2] >= first[1] + first[0]) || (first[1] >= first[2] + first[0]) || (second[0] >= second[1] + second[2]) || (second[1] >= second[2] + second[0]) || (second[2] >= second[1] + second[0]))
            answer = -1;
        else if (GeronArea(first[0], first[1], first[2]) > GeronArea(second[0], second[1], second[2]))
            answer = 1;
        else if (GeronArea(first[0], first[1], first[2]) < GeronArea(second[0], second[1], second[2]))
            answer = 2;
        else
            answer = 0;

        return answer;
        // Возвращаем большую площадь
    }
    public double GetDistance(double v, double a, double t)
    {
        return v * t + (a * t * t) / 2;
    }
    public int Task_1_3a(double v1, double a1, double v2, double a2, int time)
    {
        int answer = 0;
        if (GetDistance(v1, a1, time) == GetDistance(v2, a2, time))
        {
            answer = 0;
        }
        else if (GetDistance(v2, a2, time) < GetDistance(v1, a1, time))
        {
            answer = 1;
        }
        else
        {
            answer = 2;
        }
        return answer;
    }
    public int Task_1_3b(double v1, double a1, double v2, double a2)
    {
        int answer = 0;
        answer = 1;
        while (GetDistance(v1, a1, answer) > GetDistance(v2, a2, answer))
        {
            answer++;
        }
        return answer;
    }
    #endregion

    #region Level 2
    public void Task_2_1(int[,] A, int[,] B)
    {
        // code here
        if (A.GetLength(0) != 5 || A.GetLength(1) != 6 || B.GetLength(0) != 3 || B.GetLength(1) != 5)
        {
            return;
        }
        // create and use FindMaxIndex(matrix, out row, out column);

        FindMaxIndex(A, out int strA, out int stolbA);
        FindMaxIndex(B, out int strB, out int stolbB);
        (A[strA, stolbA], B[strB, stolbB]) = (B[strB, stolbB], A[strA, stolbA]);

    }

    public void Task_2_2(double[] A, double[] B)
    {
        // code here

        // create and use FindMaxIndex(array);
        // only 1 array has to be changed!

        // end
    }

    public void Task_2_3(ref int[,] B, ref int[,] C)
    {
        // code here
        if ((B.GetLength(0) != 5) && (C.GetLength(0) != 6) && (C.GetLength(1) != 6) && (B.GetLength(1) != 5))
        {
            return;
        }
        //  create and use method FindDiagonalMaxIndex(matrix);
        static int FindDiagonalMax(int[,] matr)
        {
            int max = -1000000;
            int index = 0;
            for (int i = 0; i < matr.GetLength(0); i++)
            {
                if (matr[i, i] > max)
                {
                    max = matr[i, i];
                    index = i;
                }
            }
            return index;
        }
        static int[,] Novmatr(int[,] matr1, int a)
        {
            int[,] matr2 = new int[matr1.GetLength(0) - 1, matr1.GetLength(1)];
            for (int i = 0; i < a; i++)
            {
                for (int j = 0; j < matr1.GetLength(1); j++)
                {
                    matr2[i, j] = matr1[i, j];
                }
            }
            for (int i = a; i < matr1.GetLength(0) - 1; i++)
            {
                for (int j = 0; j < matr1.GetLength(1); j++)
                {
                    matr2[i, j] = matr1[i + 1, j];
                }
            }
            return matr2;
        }
        B = Novmatr(B, FindDiagonalMax(B));
        C = Novmatr(C, FindDiagonalMax(C));
        // end

    }
    public void Task_2_4(int[,] A, int[,] B)
    {
        // code here

        //  create and use method FindDiagonalMaxIndex(matrix); like in Task_2_3

        // end
    }
    
    public void Task_2_5(int[,] A, int[,] B)
    {
        // code here
        int[] mas = new int[A.GetLength(1)];
        
        int a = FindMaxInColumn(A);
        int b = FindMaxInColumn(B);

        for (int j = 0; j < B.GetLength(1); j++)
        {
            mas[j] = B[b, j];
        }

        for (int j = 0; j < 6; j += 1)
        {
            B[b, j] = A[a, j];
        } 

        for (int j = 0; j < B.GetLength(1); j++)
        {
            A[a, j] = mas[j];
        }

        return;
        // end
    }

    public void Task_2_6(ref int[] A, int[] B)
    {
        // code here

        // create and use FindMax(matrix, out row, out column); like in Task_2_1
        // create and use DeleteElement(array, index);

        // end
    }
    public int CountRowPositive(int[,] matrix)
    {
        int max = -100000, c = 0, stroka = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > 0)
                    c++;
            }
            if (c > max)
            {
                max = c;
                stroka = i;
            }
            c = 0;
        }
        return c;
    }
    public int CountColumnPositive(int[,] matrix)
    {
        int max = -1000000, c = 0, stolb = 0;
        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (matrix[i, j] > 0)
                    c++;
            }
            if (c > max)
            {
                max = c;
                stolb = j;
            }
            c = 0;
        }

        return c;
    }
    public void Task_2_7(ref int[,] B, int[,] C)
    {
        // code here
        int[,] matr = new int[5, 5];
        for (int i = 0; i < CountRowPositive(B) + 1; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                matr[i, j] = B[i, j];
            }
        }
        for (int i = 0; i < 5; i++)
        {
            matr[CountRowPositive(B) + 1, i] = C[i, CountColumnPositive(C)];
        }
        for (int i = CountRowPositive(B) + 1; i < 4; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                matr[i + 1, j] = B[i, j];
            }
        }
        B = matr;
   
        // end
    }

    public void Task_2_8(int[] A, int[] B)
    {
        // code here

        // create and use SortArrayPart(array, startIndex);

        // end
    }

    public int[] Task_2_9(int[,] A, int[,] C)
    {
        
        int[] answer = new int[9];
        static int SumPositiveElementsInColumns(int[,] matrix, int a)
        {
            int sum = 0;
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (matrix[i, a] > 0)
                    sum += matrix[i, a];

            }
            return sum;
        }
        for (int i = 0; i < 4; i++)
        {
            answer[i] = SumPositiveElementsInColumns(A, i);
        }

        for (int i = 4; i < 9; i++)
        {
            answer[i] = SumPositiveElementsInColumns(C, i-4);
        }

        return answer;
    }

    public void Task_2_10(ref int[,] matrix)
    {
        // code here

        // create and use RemoveColumn(matrix, columnIndex);

        // end
    }

    public void Task_2_11(int[,] A, int[,] B)
    {
        // code here
        int a, b, a1, b1;
        FindMaxIndex(A, out a, out b);
        FindMaxIndex(B, out a1, out b1);
        int otv = A[a, b];
        A[a, b] = B[a1, b1];
        B[a1, b1] = otv;
        // end
    }
    public void Task_2_12(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxColumnIndex(matrix);

        // end
    }
    static int[,] RemoveRow(int[,] matrix, int stroka)
    {
        int[,] A = new int[matrix.GetLength(0) - 1, matrix.GetLength(1)];
        for (int i = 0; i < stroka; i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                A[i, j] = matrix[i, j];
            }
        }
        for (int i = stroka; i < matrix.GetLength(0) - 1; i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                A[i, j] = matrix[i + 1, j];
            }
        }
        return A;
    }
    public void Task_2_13(ref int[,] matrix)
    {
        int[,] A = new int[3, 5];
        int[,] B = new int[4, 5];
        static int FindMinElemStroki(int[,] A)
        {
            int min = A[0, 0];
            int stroka = 0;
            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(1); j++)
                {
                    if (A[i, j] < min)
                    {
                        min = A[i, j];
                        stroka = i;
                    }
                }
            }
            return stroka;
        }
        static int FindMaxElemStroki(int[,] A)
        {
            int max = A[0, 0];
            int stroka = 0;
            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(1); j++)
                {
                    if (A[i, j] > max)
                    {
                        max = A[i, j];
                        stroka = i;
                    }
                }
            }
            return stroka;
        }
        if (FindMaxElemStroki(matrix) != FindMinElemStroki(matrix))
        {
            B = RemoveRow(matrix, FindMaxElemStroki(matrix));
            A = RemoveRow(B, FindMinElemStroki(B));
            matrix = A;
        }
        else
        {
            B = RemoveRow(matrix, FindMaxElemStroki(matrix));
            matrix = B;
        }
    }
    public void Task_2_14(int[,] matrix)
    {
        // code here

        // create and use SortRow(matrix, rowIndex);

        // end
    }

    public int Task_2_15(int[,] A, int[,] B, int[,] C)
    {
        int answer = 0;

        // code here
        static (int, int) FindMaxElem(int[,] matrix)
        {
            int max = matrix[0, 0], stroka = 0, stolb = 0;
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    if (matrix[i, j] > max)
                    {
                        max = matrix[i, j];
                        stroka = i;
                        stolb = j;
                    }
                }
            }
            return (stroka, stolb);
        }
        static (int, int) FindMinElem(int[,] matrix)
        {
            int min = matrix[0, 0], stroka = 0, stolb = 0;
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    if (matrix[i, j] < min)
                    {
                        min = matrix[i, j];
                        stroka = i;
                        stolb = j;
                    }
                }
            }
            return (stroka, stolb);
        }
        // create and use GetAverageWithoutMinMax(matrix);
        static int GetAverageWithoutMinMax(int[,] matrix)
        {
            int a, b, a1, b1;
            (a, b) = FindMaxElem(matrix);
            (a1, b1) = FindMinElem(matrix);
            matrix[a, b] = 0;
            matrix[a1, b1] = 0;
            int c = matrix.GetLength(0) * matrix.GetLength(1) - 2;
            int sum = 0;
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    sum += matrix[i, j];
                }
            }
            return sum / c;
        }

        // 1 - increasing   0 - no sequence   -1 - decreasing

        if ((GetAverageWithoutMinMax(A) < GetAverageWithoutMinMax(B)) && (GetAverageWithoutMinMax(B) < GetAverageWithoutMinMax(C)))
        {
            answer = 1;
        }
        else if ((GetAverageWithoutMinMax(A) >= GetAverageWithoutMinMax(B)) && (GetAverageWithoutMinMax(B) >= GetAverageWithoutMinMax(C)))
            answer = -1;
        else
        {
            answer = 0;
        }
        // 1 - increasing   0 - no sequence   -1 - decreasing
        return answer;
    }

    public void Task_2_16(int[] A, int[] B)
    {
        // code here
        
        // end
    }
    public int[,] SortRowsByMaxElement(int[,] matrix1)
    {
        int[,] matrix2 = new int[matrix1.GetLength(0), matrix1.GetLength(1)];
        for (int i = 0; i < matrix2.GetLength(0); i++)
        {
            for (int j = 0; j < matrix2.GetLength(1); j++)
            {
                matrix2[i, j] = matrix1[i, j];
            }
        }
        int a = 0, b = 0;
        for (int i = 0; i < matrix1.GetLength(0); i++)
        {
            FindMaxIndex(matrix2, out a, out b);
            for (int j = 0; j < matrix1.GetLength(1); j++)
            {
                matrix1[i, j] = matrix2[a, j];
                matrix2[a, j] = -1000000;
            }
        }
        return matrix1;
    }
    public void Task_2_17(int[,] A, int[,] B)
    {
        // code here
        
        A = SortRowsByMaxElement(A);
        B = SortRowsByMaxElement(B);
    }

    public void Task_2_18(int[,] A, int[,] B)
    {
        // code here

        // create and use SortDiagonal(matrix);

        // end
    }

    public void Task_2_19(ref int[,] matrix)
    {
        //code here 
        int[,] A = new int[matrix.GetLength(0) - 1, matrix.GetLength(1)];
        int c = 0;
        for (int l = 0; l < matrix.GetLength(0); l++)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    if (matrix[i, j] == 0) c++;
                }
                if (c > 0) matrix = RemoveRow(matrix, i);
                c = 0;
            }
        }
        // end
    }
    public void Task_2_20(ref int[,] A, ref int[,] B)
    {
        // code here

        // use RemoveColumn(matrix, columnIndex); from 2_10

        // end
    }

    public void Task_2_21(int[,] A, int[,] B, out int[] answerA, out int[] answerB)
    {
        answerA = null;
        answerB = null;
        static int[] CreateArrayFromMins(int[,] matrix)
        {
            int[] mas = new int[matrix.GetLength(0)];

            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                int min = 1000000000;
                for (int j = i; j < matrix.GetLength(0); j++)
                {
                    if (min > matrix[i, j])
                    {
                        min = matrix[i, j];
                    }
                }
                mas[i] = min;
            }

            return mas;
        }
        answerA = CreateArrayFromMins(A);
        answerB = CreateArrayFromMins(B);

    }
    public void Task_2_22(int[,] matrix, out int[] rows, out int[] cols)
    {
        rows = null;
        cols = null;

        // code here

        // create and use CountNegativeInRow(matrix, rowIndex);
        // create and use FindMaxNegativePerColumn(matrix);

        // end
    }
   
    
    public void Task_2_23(double[,] A, double[,] B)
    {
        static void MatrixValuesChange(double[,] matrix)
        {

            int indx = 0;
            double[] mas = new double[matrix.GetLength(0) * matrix.GetLength(1)];

            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    mas[indx++] = matrix[i, j];
                }
            }
            for (int i = 0; i < mas.Length - 1; i++)
            {
                for (int j = 0; j < mas.Length - i - 1; j++)
                { 
                    if (mas[j] < mas[j + 1])
                    {
                        double vrem = mas[j];
                        mas[j] = mas[j + 1];
                        mas[j + 1] = vrem;
                    }
                }
            }
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    bool F = false;
                    for (int k = 0; k < 5; k++)
                    {
                        if (matrix[i, j] == mas[k])
                        {
                            F = true;
                            break;
                        }
                    }
                    if (F)
                    {
                        matrix[i, j] *= (matrix[i, j] > 0) ? 2 : 0.5;
                    }
                    else
                    {
                        matrix[i, j] *= (matrix[i, j] > 0) ? 0.5 : 2;
                    }
                }
            }
        }
        MatrixValuesChange(A);
        MatrixValuesChange(B);
    }


    public void Task_2_24(int[,] A, int[,] B)
    {
        // code here

        // use FindMaxIndex(matrix, out row, out column); like in 2_1
        // create and use SwapColumnDiagonal(matrix, columnIndex);

        // end
    }
 
    public void Task_2_25(int[,] A, int[,] B, out int maxA, out int maxB)
    {
        maxA = 0;
        maxB = 0;
        static int[] FindRowWithMaxNegativeCount(int[,] matrix)
        {
            int[] otri = new int[matrix.GetLength(0)];
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                int c = 0;
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    if (matrix[i, j] < 0)
                    {
                        c++;
                    }
                }
                otri[i] = c;
            }
            return otri;
        }
        static int otvetnegativ(int[,] matrix)
        {
            int[] negativ = FindRowWithMaxNegativeCount(matrix);
            int max = 0;
            for (int i = 0; i < negativ.Length; i++)
            {
                if (negativ[i] > negativ[max])
                {
                    max = i;
                }
            }
            return max;
        }
        maxA = otvetnegativ(A);
        maxB = otvetnegativ(B);
        // code here

        // create and use (matrix);
        // in FindRowWithMaxNegativeCount create and use CountNegativeInRow(matrix, rowIndex); like in 2_22

        // end
    }

    public void Task_2_26(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowWithMaxNegativeCount(matrix); like in 2_25
        // in FindRowWithMaxNegativeCount use CountNegativeInRow(matrix, rowIndex); from 2_22

        // end
    }
    public int strmaxindex(int[,] matrix, int row)
    {
        int ind = 0;
        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            if (matrix[row, j] > matrix[row, ind])
            {
                ind = j;
            }
        }
        return ind;
    }
    public void odd(int[,] matrix)
    {

        for (int i = 1; i < matrix.GetLength(0); i = i + 2)
        {
            int max = strmaxindex(matrix, i);
            matrix[i, max] = 0;
        }
    }
    public void even(int[,] matrix)
    {

        for (int i = 0; i < matrix.GetLength(0); i = i + 2)
        {
            int max = strmaxindex(matrix, i);
            matrix[i, max] *= (max + 1);
        }
    }
    public void Task_2_27(int[,] A, int[,] B)
    {
        // code here
        even(A);
        odd(A);
        even(B);
        odd(B);
        // create and use FindRowMaxIndex(matrix, rowIndex, out columnIndex);
        // create and use ReplaceMaxElementOdd(matrix, row, column);
        // create and use ReplaceMaxElementEven(matrix, row, column);

        // end
    }

    public void Task_2_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create and use FindSequence(array, A, B); // 1 - increasing, 0 - no sequence,  -1 - decreasing
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28b(int[] first, int[] second, ref int[,] answerFirst, ref int[,] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28c(int[] first, int[] second, ref int[] answerFirst, ref int[] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a or Task_2_28b
        // A and B - start and end indexes of elements from array for search

        // end
    }
    #endregion

    #region Level 3
    public delegate double SumFunction(int i, double x, ref int changes);
    public delegate double YFunction(double x);
   

    public void Task_3_1(ref double[,] firstSumAndY, ref double[,] secondSumAndY)
    {
        // code here
        static double sam1(int i, double x, ref int l)
        {
            if (i > 0)
            {
                l *= i;
            }
            double vrem = Math.Cos(i * x) / l;
            return vrem;
        }
        static double sam2(int i, double x, ref int k)
        {
            k *= -1;
            double vrem = k * Math.Cos(i * x) / (i * i);
            return vrem;
        }
        static double y(double x)
        {
            double y = Math.Exp(Math.Cos(x)) * Math.Cos(Math.Sin(x));
            return y;
        }
        static double y1(double x)
        {
            double y = ((x * x) - Math.Pow(Math.PI, 2) / 3) / 4;
            return y;
        }
        static double itogsam(SumFunction SumFunct, double x, int i)
        {
            double sam = 0;
            int znak = 1;
            double t = SumFunct(i, x, ref znak);
            while (Math.Abs(t) > 0.0001)
            {
                sam += t;
                t = SumFunct(++i, x, ref znak);
            }
            return sam;
        }
        static void getsam(SumFunction SumFunction, YFunction yFunction, double a, double b, double l, double[,] SamY, int one = 0)
        {
            for (int i = 0; i < (b - a) / l + 1; i++)
            {
                double x = a + i * l;
                double sum = itogsam(SumFunction, x, one);
                double y = yFunction(x);
                SamY[i, 0] = sum;
                SamY[i, 1] = y;
            }
        }
        double a = 0.1, b1 = 1, h1 = 0.1;
        double a2 = Math.PI / 5, b2 = Math.PI, h2 = Math.PI / 25;
        firstSumAndY = new double[(int)((b1 - a) / h1) + 1, 2];
        getsam(sam1, y, a, b1, h1, firstSumAndY);
        secondSumAndY = new double[(int)((b2 - a2) / h2) + 1, 2];
        getsam(sam2, y1, a2, b2, h2, secondSumAndY, 1);
        // create and use public delegate SumFunction(x) and public delegate YFunction(x)
        // create and use method GetSumAndY(sFunction, yFunction, a, b, h);
        // create and use 2 methods for both functions calculating at specific x

        // end
    }

    public void Task_3_2(int[,] matrix)
    {
        // SortRowStyle sortStyle = default(SortRowStyle); - uncomment me

        // code here

        // create and use public delegate SortRowStyle(matrix, rowIndex);
        // create and use methods SortAscending(matrix, rowIndex) and SortDescending(matrix, rowIndex)
        // change method in variable sortStyle in the loop here and use it for row sorting

        // end
    }
    public delegate void SwapDirection(double[] array);
    public double Task_3_3(double[] array)
    {
        double answer = 0;
        static void SwapLeft(double[] array)
        {
            for (int i = array.Length - 1; i > 0; i -= 2)
            {
                (array[i], array[i - 1]) = (array[i - 1], array[i]);
            }
        }
        static void SwapRight(double[] array)
        {
            for (int i = 0; i < array.Length - 1; i += 2)
            {
                (array[i], array[i + 1]) = (array[i + 1], array[i]);
            }
        }
        static double GetSum(double[] array)
        {
            double sam = 0;
            for (int i = 1; i < array.Length; i += 2)
            {
                sam += array[i];
            }
            return sam;
        }
        SwapDirection swapper = null;
        double A = array.Average();
        if (array[0] > A)
        {
            swapper = SwapRight;
        }
        else
        {
            swapper = SwapLeft;
        }
        if (swapper != null)
        {
            swapper(array);
        }
        answer = GetSum(array);




        // create and use public delegate SwapDirection(array);
        // create and use methods SwapRight(array) and SwapLeft(array)
        // create and use method GetSum(array, start, h) that sum elements with a negative indexes
        // change method in variable swapper in the if/else and than use swapper(matrix)

        // end

        return answer;
    }


    public int Task_3_4(int[,] matrix, bool isUpperTriangle)
    {
        int answer = 0;

        // code here

        // create and use public delegate GetTriangle(matrix);
        // create and use methods GetUpperTriange(array) and GetLowerTriange(array)
        // create and use GetSum(GetTriangle, matrix)

        // end

        return answer;
    }

    public void Task_3_5(out int func1, out int func2)
    {
        func1 = 0;
        func2 = 0;
        static int CountSignFlips(YFunction f, double a, double b, double h)
        {
            int znak = 0;
            double p = f(a);
            for (double x = a + h; x <= b; x += h)
            {
                double t = f(x);
                if ((p > 0 && t <= 0) || (p < 0 && t >= 0))
                {
                    znak++;
                }
                p = t;
            }
            return znak;
        }
        static double H(double x)
        {
            return x * x - Math.Sin(x);
        }
        static double H2(double x)
        {
            return Math.Exp(x) - 1;
        }
        double a1 = 0, b1 = 2, h1 = 0.1;
        double a2 = -1, b2 = 1, h2 = 0.2;
        func1 = CountSignFlips(H, a1, b1, h1);
        func2 = CountSignFlips(H2, a2, b2, h2);
        // use public delegate YFunction(x, a, b, h) from Task_3_1
        // create and use method CountSignFlips(YFunction, a, b, h);
        // create and use 2 methods for both functions

        // end
    }

    public void Task_3_6(int[,] matrix)
    {
        // code here

        // create and use public delegate FindElementDelegate(matrix);
        // use method FindDiagonalMaxIndex(matrix) like in Task_2_3;
        // create and use method FindFirstRowMaxIndex(matrix);
        // create and use method SwapColumns(matrix, FindDiagonalMaxIndex, FindFirstRowMaxIndex);

        // end
    }

    public delegate int CountPositive(int[,] matrix);
    public void Task_3_7(ref int[,] B, int[,] C)
    {
        // code here
        static int[,] InsertColumn(int[,] B, CountPositive CountstrP, int[,] C, CountPositive CountstolbP)
        {

            int indxstr = CountstrP(B);
            int indxstolb = CountstolbP(C);
            int[,] A = new int[B.GetLength(0) + 1, B.GetLength(1)];
            for (int i = 0; i < B.GetLength(1); i++)
            {
                A[indxstr + 1, i] = C[i, indxstolb];
            }
            for (int i = 0; i <= indxstr; i++)
            {
                for (int j = 0; j < B.GetLength(1); j++)
                {
                    A[i, j] = B[i, j];
                }
            }
            for (int i = indxstr + 2; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < B.GetLength(1); j++)
                {
                    A[i, j] = B[i - 1, j];
                }
            }
            return A;
        }
        B = InsertColumn(B, CountRowPositive, C, CountColumnPositive);
        // create and use public delegate CountPositive(matrix, index);
        // use CountRowPositive(matrix, rowIndex) from Task_2_7
        // use CountColumnPositive(matrix, colIndex) from Task_2_7
        // create and use method InsertColumn(matrixB, CountRow, matrixC, CountColumn);

        // end
    }
    public delegate int FindIndex(int[,] matrix);
    public void Task_3_10(ref int[,] matrix)
    {
    
    }
    public delegate int FindElementDelegate(int[,] matrix);
    public int MaxElement(int[,] matrix)
    {
        int max = int.MinValue;
        int maxindex = -1;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    maxindex = i;
                }
            }
        }
        return maxindex;
    }
    public int MinElement(int[,] matrix)
    {
        int min = int.MaxValue;
        int minindex = -1;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    minindex = i;
                }
            }
        }
        return minindex;
    }

    public void RemoveRows(ref int[,] matrix, FindElementDelegate findMax, FindElementDelegate findMin)
    {
        int max = findMax(matrix);
        int min = findMin(matrix);


        if (max == min)
        {
            matrix = RemoveRow(matrix, max);
        }
        else
        {

            matrix = RemoveRow(matrix, min);


            if (min < max)
            {
                matrix = RemoveRow(matrix, max - 1);
            }
            else
            {
                matrix = RemoveRow(matrix, max);
            }
        }
    }


    public void Task_3_13(ref int[,] matrix)
    {
        // code here
        RemoveRows(ref matrix, MaxElement, MinElement);
        // use public delegate FindElementDelegate(matrix) from Task_3_6
        // create and use method RemoveRows(matrix, findMax, findMin)

        // end
    }

    public void Task_3_22(int[,] matrix, out int[] rows, out int[] cols)
    {

        rows = null;
        cols = null;

        // code here

        // create and use public delegate GetNegativeArray(matrix);
        // use GetNegativeCountPerRow(matrix) from Task_2_22
        // use GetMaxNegativePerColumn(matrix) from Task_2_22
        // create and use method FindNegatives(matrix, searcherRows, searcherCols, out rows, out cols);

        // end
    }
    public delegate void ReplaceMaxElement(int[,] matrix);
    public void EvenOddRowsTransform(int[,] matrix, ReplaceMaxElement odd, ReplaceMaxElement even)
    {
        odd(matrix);
        even(matrix);
    }

    public void Task_3_27(int[,] A, int[,] B)
    {
        // code here
        EvenOddRowsTransform(A, odd, even);
        EvenOddRowsTransform(B, odd, even);
        // create and use public delegate ReplaceMaxElement(matrix, rowIndex, max);
        // use ReplaceMaxElementOdd(matrix) from Task_2_27
        // use ReplaceMaxElementEven(matrix) from Task_2_27
        // create and use method EvenOddRowsTransform(matrix, replaceMaxElementOdd, replaceMaxElementEven);

        // end
    }

    public void Task_3_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // create and use method FindIncreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method FindDecreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method DefineSequence(array, findIncreasingSequence, findDecreasingSequence);

        // end
    }

    public void Task_3_28c(int[] first, int[] second, ref int[] answerFirstIncrease, ref int[] answerFirstDecrease, ref int[] answerSecondIncrease, ref int[] answerSecondDecrease)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // use method FindIncreasingSequence(array, A, B); from Task_3_28a
        // use method FindDecreasingSequence(array, A, B); from Task_3_28a
        // create and use method FindLongestSequence(array, sequence);

        // end
    }
    #endregion
    #region bonus part
    public double[,] Task_4(double[,] matrix, int index)
    {
        // MatrixConverter[] mc = new MatrixConverter[4]; - uncomment me

        // code here

        // create public delegate MatrixConverter(matrix);
        // create and use method ToUpperTriangular(matrix);
        // create and use method ToLowerTriangular(matrix);
        // create and use method ToLeftDiagonal(matrix); - start from the left top angle
        // create and use method ToRightDiagonal(matrix); - start from the right bottom angle

        // end

        return matrix;
    }
    #endregion
}
