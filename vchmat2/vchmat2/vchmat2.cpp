
#include <iostream>

using namespace std;

#define eps 0.001

//вывод СЛАУ
void sysout(double** matrix, double* y, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << "(" << matrix[i][j] << "x" << j << ")";
            if (j < n - 1)
                cout << " + ";
        }
        cout << " = " << y[i] << endl;
    }
    cout << endl;
    return;
}

//диагональное преобразование
void replace(double** matrix, double* y, int n)
{
    // переставим строки так, чтобы диагоналные элементы были max
        for (int ind = 0; ind < n; ind++)
        {
            double max = matrix[0][ind]; //первый элемент текущего столбца
            int numb;
            for (int i = 1; i < n; i++)
            {
                if (matrix[i][ind] > max)
                {
                    max = matrix[i][ind];
                    numb = i;
                }
            }

            //перестановка строк,ставим на позицию ind строку, в которой ind элемент max
            if (numb != ind)
            {
                //идем по строке
                for (int i = 0; i < n; i++)
                {
                    double temp = matrix[ind][i];
                    matrix[ind][i] = matrix[numb][i];
                    matrix[numb][i] = temp;
                }

                double temp = y[ind];
                y[ind] = y[numb];
                y[numb] = temp;
            }
        }
}

//ввод СЛАУ
void sysin(double** matrix, double* y, int n)
{
    cout << "Введите коэффициенты при неизвестных:" << endl;;
    for (int i = 0; i < n; i++)
    {
        matrix[i] = new double[n];
        for (int j = 0; j < n; j++)
            cin >> matrix[i][j];
    }
    cout << endl;

    cout << "Введите коэффициенты правой части уравнений:" << endl;;
    for (int i = 0; i < n; i++)
        cin >> y[i];
    cout << endl;
    return;
}

//метод Гаусса без выбора главного элемента
double* gauss(double** matrix, double* y, int n)
{
    double* x;
    x = new double[n];

    // переставим строки так, чтобы диагоналные элементы были не 0
    for (int ind = 0; ind < n; ind++)
    {
        double max = matrix[0][ind]; //первый элемент текущего столбца
        int numb;
        for (int i = 1; i < n; i++)
        {
            if (matrix[i][ind] != 0)
            {
                //max = matrix[i][ind];
                numb = i;
            }
        }

        //перестановка строк,ставим на позицию ind строку, в которой ind элемент max
        if (numb != ind)
        {
            //идем по строке
            for (int i = 0; i < n; i++)
            {
                double temp = matrix[ind][i];
                matrix[ind][i] = matrix[numb][i];
                matrix[numb][i] = temp;
            }

            double temp = y[ind];
            y[ind] = y[numb];
            y[numb] = temp;
        }
    }

    //идем по переменным
    for (int ind = 0; ind < n; ind++)
    {
        //приведем расширенную матрицу к ступенчатому виду
        //идем по строкам,начиная со следующей после выбранной
        for (int i = ind + 1; i < n; i++)
        {
            double mult = -matrix[i][ind] / matrix[ind][ind];
            if (matrix[i][ind] != 0)
            {
                for (int j = ind; j < n; j++)
                    matrix[i][j] += matrix[ind][j] * mult;
            }
            else
                continue;
            y[i] += y[ind] * mult;
        }
    }

    cout << "Ступенчатая матрица:" << endl;
    sysout(matrix, y, n);

    //обратная подстановка
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = y[i] / matrix[i][i];
        for (int j = 0; j < i; j++)
            y[j] = y[j] - matrix[j][i] * x[i];
    }

    cout << "Решение системы:" << endl;
    return x;
}

//критерий окончания итераций
bool end(double* x, double* xPrev, int n, double** matrix)
{
    double am, max, norm = 0;

    //ищем m-норму матрицы
    for (int i = 0; i < n; i++)
    {
        am=0, max = 0;
        for (int j = 0; j < n; j++)
            am += abs(matrix[i][j]);
        if (am > max)
            max = am;
    }
    am = max;
    cout << "am = " << am;
    for (int i = 1; i < n; i++) 
        norm += pow((x[i] - xPrev[i]), 2);
    norm = sqrt(norm);
    return (norm <= (1 - am) * eps / am);
}

//вывод приближений
void print(double* x, int n)
{
    for (int i = 0; i < n; i++)
        cout <<"x"<<i<<":   "<< " = " << x[i] << endl;
    return;
}

//метод Зейделя
void zeidel(double** matrix, double* y, int n)
{
    double* x0;  //массив начальных приближений k
    double* b;  //свободный член в правой части экв.системы
    double** matrix2;  //коэф.при х в экв.системе

    double* x0Next;  //массив начальных приближений k+1
    int k = 1;  //количество итераций

    double am, al, max1, max2;

    x0 = new double[n];
    b = new double[n]; 
    matrix2 = new double* [n]; 
    x0Next = new double[n];
   
    //делаем перестановку, чтобы диагональные элементы были max
    replace(matrix, y, n);
    sysout(matrix, y, n);

    //приводим систему к эквивалентному виду
    for (int i = 0; i < n; i++)
    {
        matrix2[i] = new double[n];
        cout << "x" << i << " = "; 
        b[i] = y[i] / matrix[i][i];  
        cout << b[i]<<" + ";

        for (int j = 0; j < n; j++) 
        {
            if (j != i) 
                matrix2[i][j] = -matrix[i][j] / matrix[i][i];
            else
                matrix2[i][j]=0;
            cout << "(" << matrix2[i][j] << "x" << j << ")";
            if (j < n - 1)
                cout << " + ";
        }
        cout << endl;
    }

    //критерий сходимости
    //ищем m-норму матрицы
    for (int i = 0; i < n; i++)
    {
        am = 0, max1 = 0;
        for (int j = 0; j < n; j++)
            am += abs(matrix2[i][j]);
        if (am > max1)
            max1 = am;
    }
    am = max1;
    //ищем l-норму матрицы
    for (int i = 0; i < n; i++)
    {
        al = 0, max2 = 0;
        for (int j = 0; j < n; j++)
            al += abs(matrix2[j][i]);
        if (al > max2)
            max2 = al;
    }
    al = max2;
    if ((am < 1) || (al < 1))
        cout << "Процесс итерации сходится.";
    else
    {
        cout << "Процесс итерации не сходится.";
        return;
    }

   cout << endl << "Началное приближение x0" << endl;
   for (int i = 0; i < n; i++)
   {
       x0Next[i] = b[i];
       cout <<"x"<<i<<" = "<< x0Next[i] << endl;
   }

   do
    {
       for (int i = 0; i < n; i++)
           x0[i] = x0Next[i];
       
       //по строкам для каждой переменной ищем (k+1)-ое приближение
        for (int i = 0; i < n; i++)
        {
            x0Next[i] = b[i];
            int j;
            double sum1 = 0, sum2 = 0;

            for (j = 0; j < i; j++)
                sum1 += matrix2[i][j] * x0Next[j];

            for (int ii = j; ii < n; ii++)
                sum2 += matrix2[i][ii] * x0[ii];

            x0Next[i] += sum1 + sum2;
        }
        cout <<endl<< "Итерация " << k << ":" << endl;
        print(x0Next, n);
        k++;
    } while (!end(x0Next,x0,n,matrix2));

    return;
}

int main()
{
    setlocale(LC_ALL, "Russian");

    double** matrix;
    double* y;
    int n;
    double* x;

    cout << "Введите количество уравнений: ";
    cin >> n;
    cout << endl;

    matrix = new double* [n];
    y = new double[n];

    sysin(matrix, y, n);
    cout << "Система уравнений: " << endl;
    sysout(matrix, y, n);

    //Решение методом Гаусса
    x = gauss(matrix, y, n);
    for (int i = 0; i < n; i++)
        cout << "x" << i << "= " << x[i] << endl;
  
    //решение методом Зейделя
    //zeidel(matrix, y, n);
}

