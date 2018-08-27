#pragma once
#include <iostream>
using namespace std;

struct Point
{
    float coords[2];
    int index;
};

struct Comparator
{
    static int coord;
    static int compare_point(const void *a, const void *b);
};

bool check_sorted(Point *points, int n);

float generate_float_number();

void fill_Point_array(Point *point_array, long number_1, long number_2, long round_total);

void print_points(Point *points, long n);