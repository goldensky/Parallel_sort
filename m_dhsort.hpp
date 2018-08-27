
#include "Point.h"

void down_heap(Point *a, long k, long n)
{
    Point elem;
    long child;
    elem = a[k];
    
    while (k <= n / 2)
    {
        child = 2 * k;
        
        if (child < n && Comparator::compare_point(&a[child], &a[child + 1]) < 0)
            child++;
        if (Comparator::compare_point(&elem, &a[child]) >= 0)
            break;
        
        a[k] = a[child];
        k = child;
    }
    a[k] = elem;
}

void heap_sort(Point *a, long size)
{
    long i;
    Point temp;
    
    for (i = size - 1; i >= 0; --i)
        down_heap(a, i, size - 1);
    
    for (i = size - 1; i > 0; --i)
    {
        temp = a[i];
        a[i] = a[0];
        a[0] = temp;
        
        down_heap(a, 0, i - 1);
    }
}

void merge(Point *a, Point *b, long size1, long size2)
{
    long pos1 = 0, pos2 = 0, pos3 = 0;
    
    Point *temp = new Point[size1 + size2];
    
    while (pos1 < size1 && pos2 < size2)
    {
        if (Comparator::compare_point(&a[pos1], &b[pos2]) < 0)
            temp[pos3++] = a[pos1++];
        else
            temp[pos3++] = b[pos2++];
    }
    
    while (pos2 < size2)
        temp[pos3++] = b[pos2++];
    while (pos1 < size1)
        temp[pos3++] = a[pos1++];
    
    for (pos3 = 0; pos3 < size1; ++pos3)
        a[pos3] = temp[pos3];
    
    pos2 = 0;
    for (pos3 = size1; pos3 < size1 + size2; ++pos3)
        b[pos2++] = temp[pos3];
    
    delete[] temp;
}

void dhsort(Point *a, long size)
{
    if (size == 1)
        return;
    
    if (size == 2)
    {
        if (Comparator::compare_point(&a[0], &a[1]) > 0)
        {
            Point tmp = a[0];
            a[0] = a[1];
            a[1] = tmp;
        }
        return;
    }
    
    long half = size / 2, rest = size - size / 2;
    if (size > 100000)
    {
        dhsort(a, half);
        dhsort(&a[half], rest);
        merge(a, &a[half], half, rest);
    }
    else
        heap_sort(a, size);
}

void dhsort_par(Point *a, long size, int threads)
{
    if (size == 1)
        return;
    
    if (threads == 1)
    {
        dhsort(a, size);
        return;
    }
    
    long half = size / 2, rest = size - size / 2;
#pragma omp parallel sections
    {	
#pragma omp section	
        dhsort_par(a, half, threads / 2);
#pragma omp section
        dhsort_par(&a[half], rest, threads - threads / 2);
    }
    
    merge(a, &a[half], half, rest);
    
    return;
}
