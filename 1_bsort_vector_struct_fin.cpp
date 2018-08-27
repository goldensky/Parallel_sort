#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>

const int TEST_NUMBER_MAX = 24;  // 24

using namespace std;

//установка флага в true  показывает вывод на печать промежуточных результатов
int DEBUG =  false; //true;
int DEBUG1 = false;
int DEBUG2 = false;

struct Comparator
{
    int a;
    int b;
};

// печать результата (по заданному образцу)
void print_result(vector<Comparator> Comparator_array_vector,  int numbers, int n_comp,  int n_tact)
{
    cout << numbers << ' ' << 0 << ' ' << 0 << endl;
    
    if (numbers == 1)
    {
        cout << 0 << endl << 0 << endl;
    }
    else if (numbers == 2)
    {
        cout << " " << 0 << " " << 1 << endl;
        cout << 1 << endl << 1 << endl;
    }
    
    else
    {
        for(vector<Comparator>:: iterator it = Comparator_array_vector.begin(); it != Comparator_array_vector.end(); ++it)
        {
            cout << " " << it -> a << " " << it -> b << endl;
        }
        cout << n_comp << endl;
        cout << n_tact << endl;
    }
}//print_result

// операция сравнения - обмена для сортировки Бэтчера
void compexch(vector<int> &A, int elem1 ,  int elem2)
{
    if(A[elem1] > A[elem2])
    {
        int temp = A[elem1];
        A[elem1] = A[elem2];
        A[elem2] = temp;
    }
}//compexch

void print_array(vector<int> A)  // not need
{
    for(size_t i = 0; i < A.size(); i++)
    {
        cout << A[i] << ' ' ;
    }
    cout << endl;
}

int max(int a, int b)
{
    return (a > b) ? a : b;
}

int find_tacts(vector<int> tact_array)
{
    int m = tact_array[1];
    for(size_t i = 2; i < tact_array.size(); ++i)
    {
        if (m < tact_array[i])
        {
            m = tact_array[i];
        }
    }
    
    
    return m;
    
}

// Сортировка Бэтчера из Искусства Программирования том 3 стр 132
void sort(vector<int> &A, int numbers,  bool test_print = true)
{
    int n_comp = 0;
    int n_tact = 0;
    
    vector<Comparator> Comparator_array_vector;
    vector<int> tact_array(numbers);
        
    int t = ceil(log2(numbers));
    int p = pow(2, t - 1);
    
    for(;p > 0; p = (int)floor(p / 2))
    {
        int r = 0;
        int d = p;
        int q = pow(2, t - 1);
        
        while(q >= p)
        {
            for(int i = 0; i < numbers - d; ++i)
            {
                if((i & p) == r)
                {
                    Comparator comp;
                    comp.a = i;
                    comp.b = i + d;
                    Comparator_array_vector.push_back(comp);
                    
                    
                    if(A.size() != 0)
                    {
                        compexch(A, i, i+d);
                    }
                    
                    n_comp++;
                }
            }
            
            d = q - p;
            q = q / 2;
            r = p;
        }
        
    }
    
    
    //Нахождение количества тактов работы сортировки Бэтчера
    ///////
   
    for(int i = 0; i < numbers; i++)
    {
        tact_array[i] = 0;
    }
   
    for(vector<Comparator>:: iterator it = Comparator_array_vector.begin(); it != Comparator_array_vector.end(); ++it)
    {
        int max_value = max(tact_array[it->a], tact_array[it->b]) + 1;
        tact_array[it->a] = tact_array[it->b] = max_value;
    }    
 
    n_tact = find_tacts(tact_array); // !
   
    if (test_print)
    {
        print_result(Comparator_array_vector, numbers, n_comp, n_tact);
    }
    
}//sort

// проверка отсортированности полученного массива
bool is_sorted(vector<int> &A)
{
    bool OK = true;
    for (size_t i = 0; i < A.size() - 1; ++i)
    {
        if(A[i] > A[i+1])
        {
            cout << "NONONONONONON" << endl;
            OK = false;
            break;
        }
        
    }
    
    if(DEBUG)
    {
        // проверочная печать отсортированного массива
        cout << "sorted array " << endl;
        for (size_t i = 0; i < A.size(); ++i)
        {
            cout << A[i] << ' '  ;
        }
        cout << endl << "--------------" << endl;
    }
        
    return OK;
}//is_sorted


// создает массив из нулей и единиц  -  все перестановки -  их  2**numbers
// затем этот массив сортируем для проверки работоспособности программы
void create_array(vector<int> &A, unsigned mask)
{
    for(size_t i = 0; i < A.size(); ++i)
    {
        A[i] = mask & 1;
        mask >>= 1;
    }
    
    if(DEBUG)
    {
        cout << endl << "create array" << endl;
        for(size_t i = 0; i < A.size(); ++i)
        {
            cout << A[i] << ' ';
        }
        
        cout << "\n---------------" << endl << endl;
    }
    
    
}// create_array

// вспомогательная функция для проверки работоспособности программы
// используется в  test_n()
bool check_function(int numbers)
{
    vector<int> A(numbers);
    int top = 1 << numbers;
    bool OK = true;
    
    for (int mask = 0; mask < top; mask++)
    {
        create_array(A, mask);
        //count_all_permutations++;
        sort(A, numbers,  false);
        
        if (!is_sorted(A))
        {
            cout << "not_sorted " << endl;
            OK = false;
            return false;
        }
    }
    
    if(DEBUG1)
        cout << "OK = " << OK << endl;
    
    return OK;
    
}//check_function

// итоговая тестовая функция для проверки  на 24 элементах
// можно произвольное количество  - TEST_NUMBER записан в начале программы
void test_n(const int n)
{
    for (int numbers = 1; numbers <= n; ++numbers)
    {
		cout << "Check for length: "<<numbers<<" ..."<<endl;
        if(!check_function(numbers))
        {
            cerr << "Wrong result for " << numbers << endl;
            break;
        } else {
			cout << "Ok."<< endl;
		}
        if(DEBUG)
            cout << endl << "-------------------------------------------------" << endl;
    }
}// test_24

int main(int argc, char *argv[])
{
    
    int numbers = 0;
    if (argc > 1)
    {
        numbers = atoi(argv[1]);
    }
    else
    {
        cout << "Please input number of elements" << endl;
        return 0;
    }
    
    
    if (numbers > 10000)
    {
        cout <<   "Too many elements. Please input number of elements less then 10000" << endl;
        return 0;
    }
    
    
    vector<int> A;
    

    
    if (DEBUG2)
    {
        printf("+++++++++++++++++++++++++++++++++\n");
        int OK = check_function(numbers);
        printf("OK = %d \n", OK);
    }
    
    unsigned long long int start_time = clock();
    
    sort(A, numbers);
    
    unsigned long long int end_time = clock();
    cout <<  endl << "time sort = " << end_time - start_time << " for numbers =  " << numbers << endl;
    
    
    //задана проверка для перестановок нулей и единиц длины 24 в соответствии с принципом нулей и единиц
    //если неадаптивная программа выдает отсортированный выход, когда входы состоят только из 0 и 1,
    //то она делает то же самое, когда входами являются произвольные ключи
    
    start_time = clock();
    
    
    if (argc == 2)
    {
        printf("The defaul value for checking  n = %d\n", TEST_NUMBER_MAX);
        test_n(TEST_NUMBER_MAX);
    }
    
    
    
    if (argc > 2 ) {
		const int check_till = atoi(argv[2]);
		if (check_till < 0 || check_till > TEST_NUMBER_MAX)
        {
			printf("Expected number between 1 and %d. The defaul value for checking  n = %d\n", TEST_NUMBER_MAX, TEST_NUMBER_MAX);
            test_n(TEST_NUMBER_MAX);
		}
        else
        {
			test_n(check_till);
		}
	}
    
    end_time = clock();
    cout <<  endl << "time test_n  = " << end_time - start_time << endl;
    
    
    
    return 0;
}
