#include <iostream>
#include <time.h>

#include "m_dhsort.hpp"

using namespace std;



//    mpicxx  -Wall -o task2   posl_new.cpp  Point.cpp -fopenmp -pthread
//    mpirun --allow-run-as-root  -n 1 ./task2 1000 1000 1



int main(int argc, char **argv)
{
    int coordinate_for_compare = 0;
    long number_1 = 0, number_2 = 0;
    coordinate_for_compare = 0;
    
    if(argc < 3)
    {
        cout <<  "Please input at least  dimension of network: number_1   number_2 " << endl;
        cout <<  "or input dimension of network  and coord number (0 or 1): " ;
        cout << "  number_1   number_2   coord" << endl;
        return 0;
    }
    
    else if(argc == 3)
    {
        number_1 = atoi(argv[1]);
        number_2 = atoi(argv[2]);
        if (number_1 <= 0)  return 0;
        if (number_2 <= 0)  return 0;
        coordinate_for_compare = 0;
    }
    
    else if (argc >= 4)
    {
        number_1 = atoi(argv[1]);
        number_2 = atoi(argv[2]);
        if (number_1 <= 0)  return 0;
        if (number_2 <= 0)  return 0;
        
        
        coordinate_for_compare = atoi(argv[3]);
        if (coordinate_for_compare != 0 && coordinate_for_compare != 1)
        {
            coordinate_for_compare = 0;
        }
    }

    Comparator::coord = coordinate_for_compare;
    
    Point *points = new Point[number_1 * number_2];
    
    fill_Point_array(points, number_1, number_2,  number_1 * number_2);

    
    clock_t sort_time = clock();
    
    //qsort(points, number_1 * number_2, sizeof(*points), Comparator::compare_point);
    
    //heap_sort(points, number_1 * number_2);
    
    dhsort(points, number_1 * number_2);
    
    sort_time = clock() - sort_time;
    
    
    
    //print_points(points, number_1 * number_2);
    
    bool sorted = check_sorted(points, number_1 * number_2);
    
    cout << "number_1 = " << number_1 << endl;
    
    cout << "number_2 = " << number_2 << endl;
    cout << "sorted = " << sorted << endl;
    
    cout << "sizeof(Point) = " << sizeof(Point) << endl;
    
    cout << "Sort time (sec): " << sort_time * 1.0 / CLOCKS_PER_SEC << endl;
    
    delete[] points;
}
