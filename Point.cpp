#include "Point.h"
#include <cmath>

int Comparator::coord = 0;

int Comparator::compare_point(const void *ap, const void *bp)
{
    Point *a = (Point*)ap;
    Point *b = (Point*)bp;
    return a->coords[coord] < b->coords[coord] ? -1 : a->coords[coord] == b->coords[coord] ? 0 : 1;
    
}//



bool check_sorted(Point *points, int n)
{
    for (int i = 0; i < n - 1; i++)
        if (Comparator::compare_point(&points[i+1], &points[i]) < 0 )
            return false;
    return true;
}



//
float generate_float_number()
{
    return (float)  (rand() / (RAND_MAX / 1000000.0));
}//generate_float_number



void fill_Point_array(Point *point_array, long number_1, long number_2, long round_total)
{
    long i, j;
    long index;
    
    for (i = 0; i < number_1; i++)
    {
        for (j = 0; j < number_2; j++)
        {
            index = i * number_2 + j;
            point_array[index].coords[0] = generate_float_number() * i + j;
            point_array[index].coords[1] = generate_float_number() * j + i;
            
            point_array[index].index = index;
            //cout << "index = " << index << endl;
            //cout << "  0      i  j         "  <<  i << ' ' << j << ' ' << index << endl;
        }
    }
    
    // false point
    index++;
    for (; index < round_total; index++)
    {
        point_array[index].index = -1;
        point_array[index].coords[0] = -1;
        point_array[index].coords[1] = -1;
        //cout << "  1      i  j             "  <<  i << ' ' << j << ' ' << index << endl;
    }
    
}//  fill_Point_array




//
void print_points(Point *points, long n)
{
    for (long i = 0; i < n; i++)
    {
        cout << points[i].index << ' ';
        cout << points[i].coords[0] << ' ' << points[i].coords[1] << ' ' << endl;
    }
}// print_points

