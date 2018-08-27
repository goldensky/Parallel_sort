#include <omp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "mpi.h"
#include "m_dhsort.hpp"
#include "Point.h"
#include <float.h>

using namespace std;


//  mpisubmit.bg -n 128 --w 00:01:00 -e "OMP_NUM_THREADS=4" --stdout 500x500_128.txt task2 500 500 0

//  
//  mpixlcxx_r -qsmp=omp -qarch=450d -qtune=450 -o prog_name

//  mpisubmit.bg --nproc 128 --mode smp  ./task2bg.exe [--your program args]

mpisubmit.bg --nproc 128 --mode smp  ./taskD
mpisubmit.bg --nproc 4 --mode smp  ./taskD


//    mpicxx  -Wall -o    task2par   parallel2.cpp  Point.cpp -fopenmp -pthread
//    mpirun --allow-run-as-root  -n 2 ./task2par    1000 1000 1




mpisubmit.bg -n 128  -e "OMP_NUM_THREADS=4"  taskD 100 100



mpicxx  -Wall -o taskD   task2.cpp   -fopenmp -pthread



vector<pair<int, int> > comparators;
vector<int> lines;
MPI_Datatype MPI_POINT_TYPE;


void swap_ptrs(void *ptr1, void *ptr2)
{
    void **pptr1 = (void **)ptr1;
    void **pptr2 = (void **)ptr2;
    void *temp = *pptr1;
    *pptr1 = *pptr2;
    *pptr2 = temp;
    
}

void batcher_merge(const vector<int> &n1, const vector<int> &n2)
{
    if (n1.size() == 0 || n2.size() == 0)
        return;
    
    if (n1.size() + n2.size() == 2)
    {
        comparators.push_back(make_pair(n1[0], n2[0]));
        lines[n1[0]] = lines[n2[0]] = max(lines[n1[0]], lines[n2[0]]) + 1;
        return;
    }
    
    vector<int> n1_odd, n2_odd, n1_even, n2_even;
    
    for (int i = 0; i < n1.size(); ++i)
    {
        if (i % 2 == 0)
            n1_even.push_back(n1[i]);
        else
            n1_odd.push_back(n1[i]);
    }
    
    for (int i = 0; i < n2.size(); ++i)
    {
        if (i % 2 == 0)
            n2_even.push_back(n2[i]);
        else
            n2_odd.push_back(n2[i]);
    }
    
    batcher_merge(n1_odd, n2_odd);
    batcher_merge(n1_even, n2_even);
    
    vector<int> res(n1.begin(), n1.end());
    res.insert(res.end(), n2.begin(), n2.end());
    
#pragma omp parallel
    {
        vector<pair<int, int> > temp;
#pragma omp for nowait schedule(static)
        for (int i = 1; i < res.size() - 1; i += 2)
        {
            temp.push_back(make_pair(res[i], res[i + 1]));
            lines[res[i]] = lines[res[i + 1]] = max(lines[res[i]], lines[res[i + 1]]) + 1;
        }
        
#pragma omp for schedule(static) ordered
        for (int i = 0; i < omp_get_num_threads(); ++i)
        {
#pragma omp ordered
            comparators.insert(comparators.end(), temp.begin(), temp.end());
        }
    }
}

void batcher(const vector<int> &n)
{
    if (n.size() == 1)
        return;
    
    vector<int> n1(n.begin(), n.begin() + n.size() / 2);
    vector<int> n2(n.begin() + n.size() / 2, n.end());
    
    batcher(n1);
    batcher(n2);
    batcher_merge(n1, n2);
}


bool check_sorted_total(Point *points, int items_per_processor, int rank, int processors, MPI_Datatype MPI_POINT)
{
    bool local_sorted = check_sorted(points, items_per_processor);
    bool borders_sorted = true;
    
    int right_rank = (rank + 1) % processors;
    int left_rank = (rank - 1) % processors;
    
    Point *left_max = new Point();
    MPI_Sendrecv(
                 points + items_per_processor - 1, 1, MPI_POINT, right_rank, 1,
                 left_max, 1, MPI_POINT, left_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE
                 );

    if (rank > 0 && Comparator::compare_point(&points[0], *&left_max) < 0)
        
    {
        borders_sorted = false;
    }
    int local_borders_sorted = local_sorted && borders_sorted;
    
    cout << " rank = " << rank << "  " << "   local_sorted = " << local_sorted << "  borders_sorted " << borders_sorted;
    cout <<  "    local_borders_sorted  =   " << local_borders_sorted    <<   endl;
    
    
    int global_sorted;
    MPI_Allreduce(&local_borders_sorted, &global_sorted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    return global_sorted == processors;
}


void create_datatype()
{
    Point p;
    
    // new   MPI  type
    
    const int n = 2;
    
    int blocks[n] = {2, 1};
    MPI_Aint offset[n];
    
    MPI_Datatype data_type[n] = {MPI_FLOAT, MPI_INT};
    
    offset[0] = offsetof(Point, coords);
    offset[1] = offsetof(Point, index);
    
    MPI_Type_create_struct(n, blocks, offset, data_type, &MPI_POINT_TYPE);
    MPI_Type_commit(&MPI_POINT_TYPE);
    
}


Point create_false_point()
{
    Point fp;
    fp.coords[0] = FLT_MAX;
    fp.coords[1] = FLT_MAX;
    fp.index = -1;
    return fp;
}


void generate_points(Point *points, long n, int rank, long total_items)
{
    for(int j = 0; j < n; j++)
    {
        int index = n * rank + j;
        if (index >= total_items)
        {
            points[j] = create_false_point();
        }
        else {
            points[j].coords[0] = generate_float_number();
            points[j].coords[1] = generate_float_number();
            points[j].index = index;
        }
    }
}


int main(int argc, char **argv)
{
    int coordinate_for_compare = 0;
    long number_1 = 0, number_2 = 0;
    int num_threads = 1;

    if(argc < 3)
    {
        cout <<  "Please input at least  dimension of network: number_1   number_2 " << endl;
        cout <<  "or input dimension of network, threads number and coord number (0 or 1): " ;
        cout << "  number_1   number_2   threads   crd" << endl;
        return 0;
    }
    
    else if(argc == 3)
    {
        number_1 = atoi(argv[1]);
        number_2 = atoi(argv[2]);
        num_threads = 1;
        coordinate_for_compare = 0;
        
    }
    
    else if (argc == 4)
    {
        number_1 = atoi(argv[1]);
        number_2 = atoi(argv[2]);
        num_threads = atoi(argv[3]);
        if(num_threads <= 0)
            num_threads = 1;
        coordinate_for_compare = 0;
    }
    
    
    else if(argc >= 5)
    {
        number_1 = atoi(argv[1]);
        number_2 = atoi(argv[2]);
        num_threads = atoi(argv[3]);
        if(num_threads <= 0)
            num_threads = 1;
        
        coordinate_for_compare = atoi(argv[4]);
        if (coordinate_for_compare != 0 && coordinate_for_compare != 1)
        {
            coordinate_for_compare = 0;
        }
    }
    
    if (number_1 <= 0)  return 0;
    if (number_2 <= 0)  return 0;
    
    Comparator::coord = coordinate_for_compare;
    
    int rank;
    int processors;
    
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    
    MPI_Comm_size(MPI_COMM_WORLD, &processors);
    MPI_Comm_rank(comm, &rank);
    
    if (processors < 1)
    {
        MPI_Finalize();
        return 0;
    }
    
    
    omp_set_num_threads(num_threads);
    num_threads = omp_get_max_threads();

    vector<int> array;
    for (int i = 0; i < processors; ++i)
        array.push_back(i);
    lines = vector<int>(processors, 0);
    batcher(array);
    
    /*
     for (pair<int, int> elem : comparators)
     cout << elem. first << ' ' << elem. second << '\n';
     */
    
    MPI_Barrier(comm);
    
    long items_per_processor = ceil(number_1 * number_2 / (float)processors);
    Point *points = new Point[items_per_processor];

    generate_points(points, items_per_processor, rank, number_1 * number_2) ;

    MPI_Status status;
    create_datatype();
    
    
    Point *tmp1 = new Point[items_per_processor];
    Point *tmp2 = new Point[items_per_processor];

    
    double sort_time = MPI_Wtime(), total = sort_time;
    
    
    dhsort_par(points, items_per_processor, num_threads);
    
    sort_time = MPI_Wtime() - sort_time;
    
    double btime = MPI_Wtime();
    
    for (int i = 0; i < comparators.size(); ++i)
    {
        pair<int, int> comp = comparators[i];
        if (rank == comp.first)
        {
            MPI_Send(points, items_per_processor, MPI_POINT_TYPE, comp.second, 0, MPI_COMM_WORLD);
            MPI_Recv(tmp2, items_per_processor, MPI_POINT_TYPE, comp.second, 0, MPI_COMM_WORLD, &status);
            
            long ind_rez = 0;
            long ind_cur = 0;
            for (long ind_tmp = 0; ind_tmp < items_per_processor; ++ind_tmp)
            {
                Point res = points[ind_rez];
                Point cur = tmp2[ind_cur];
                if (Comparator::compare_point(&res, &cur) < 0)
                {
                    tmp1[ind_tmp] = res;
                    ind_rez++;
                }
                else {
                    tmp1[ind_tmp] = cur;
                    ind_cur++;
                }
            }
            swap_ptrs(&points, &tmp1);
        }
        
        if (rank == comp.second)
        {
            MPI_Recv(tmp2, items_per_processor, MPI_POINT_TYPE, comp.first, 0, MPI_COMM_WORLD, &status);
            MPI_Send(points, items_per_processor, MPI_POINT_TYPE, comp.first, 0, MPI_COMM_WORLD);
            
            long ind_rez = items_per_processor - 1;
            long ind_cur = items_per_processor - 1;
            for (long ind_tmp = items_per_processor - 1; ind_tmp >= 0; --ind_tmp)
            {
                Point res = points[ind_rez];
                Point cur = tmp2[ind_cur];
                if (Comparator::compare_point(&res, &cur) > 0)
                {
                    tmp1[ind_tmp] = res;
                    ind_rez--;
                }
                else
                {
                    tmp1[ind_tmp] = cur;
                    ind_cur--;
                }
            }
            swap_ptrs(&points, &tmp1);
        }
    }
    
    double temp = MPI_Wtime();
    btime = temp - btime;
    total = temp - total;
    
    double max_total, sort_max, btime_max;
    
    MPI_Reduce(&total, &max_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sort_time, &sort_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&btime, &btime_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    bool sorted = check_sorted_total(points, items_per_processor, rank, processors, MPI_POINT_TYPE);
    
    cout <<  "! sorted = " <<  sorted  << endl;
    
    MPI_Barrier(comm);
    
    if (rank == 0)
    {
        cout << endl << endl << "Print Report" << endl;
        cout << endl << endl << "items_per_processor" << items_per_processor << endl;
        cout <<  "number_1 :   " << number_1 << endl;
        cout <<  "number_2  :  " << number_2 << endl;
        cout << "Number of elements: " << number_1 * number_2 << endl;
        cout << "coordinate_for_compare: " << coordinate_for_compare << endl;
        cout << "Processors: " << processors << endl;
        cout << "Threads: " << num_threads << endl;
        cout << "Batcher steps: " << *max_element(lines.begin(), lines.end()) << endl;
        cout << "Comparators count: " << comparators.size() << endl;
        cout << "Total time (sec): " << max_total << endl;
        cout << "Sort time (sec): " << sort_max << endl;
        cout << "Batcher time (sec): " << btime_max << endl;
    }
    
    
    delete[] tmp1;
    delete[] tmp2;
    delete[] points;
    
    MPI_Type_free(&MPI_POINT_TYPE);
    MPI_Finalize();
    
    return 0;
}