#include "omp.h"
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

typedef long long data_t;
#define N 10000000
typedef long long index_t;
typedef bool (*predicate)(data_t);
typedef struct filter_ret_t
{
    data_t* filtered_array;
    index_t filtered_array_len;
}filter_ret_t;

data_t *prefix_sum(data_t *array, index_t n)
{
    data_t *x = malloc(sizeof(data_t) * n);
    data_t *t = malloc(sizeof(data_t) * n);
    index_t i,j;
    #pragma omp parallel for private(i)
    for (i = 0; i < n; i++)
        x[i] = array[i];
    for (j = 0; j < log2(n); j++)
    {
    #pragma omp parallel private(i)
        {
            #pragma omp for
            for (i = 1 << j; i < n; i++)
                t[i] = x[i] + x[i - (1 << j)];
            #pragma omp for
            for (i = 1 << j; i < n; i++)
                x[i] = t[i];
        }
    }
    free(t);
    return x;
}

filter_ret_t filter(data_t *array, index_t length, predicate p)
{
    data_t *filtered, *bitmap, *bitsum;
    bitmap = malloc(sizeof(data_t) * length);
    #pragma omp parralel for
    for(int i=0;i<length;i++)
    {
        bitmap[i] = p(array[i]);
    }
    bitsum = prefix_sum(bitmap, length);
    index_t filtered_length = bitsum[length - 1];
    filtered = malloc(sizeof(data_t) * filtered_length);
    if(bitsum[0] > 0) // edge index
    {
        filtered[0] = array[0];
    }
    #pragma omp parallel for
    for (int i = 1; i < length; i++)
    {
        if(bitsum[i] > bitsum[i-1])
        {
            filtered[bitsum[i]-1] = array[i];
        }
    }
    free(bitsum);
    free(bitmap);
    filter_ret_t ret = {filtered, filtered_length};
    return ret;
}

bool even(data_t a)
{
    return a % 2 == 0;
}

void test_psum()
{
    bool assert_psum = true;
    data_t *a = malloc(sizeof(data_t) * N);
    for (long i = 0; i < N; i++)
        a[i] = i + 1;
    data_t *a_psum = prefix_sum(a, N);
#pragma omp parallel for reduction(&: assert_psum)
    for (index_t i = 1; i <= N; i++)
    {
        data_t expected = (i * (i + 1)) / 2;
        bool current_test = (a_psum[i - 1] == expected);
        if (!current_test)
        {
            printf("test %d failed. expected %d, actual %d\n", i, expected, a_psum[i]);
        }
        assert_psum = assert_psum && current_test;
    }
    printf("assert_psum=%s\n", assert_psum ? "True" : "False");
    free(a_psum);
    free(a);
}

// void test_filter()
// {
//     data_t *a = malloc(sizeof(data_t) * N);
//     for (long  i = 0; i < N; i++)
//         a[i] = i+1;
//     data_t *a_psum = prefix_sum(a, N);
//     filter_ret_t a_filtered = filter(a, N, even);
//     bool assert_filtered = true;
//     #pragma omp parallel for reduction(&: a_filtered)
//     for (index_t i = 1; i <= a_filtered.filtered_array_len; i++)
//     {
//         data_t expected = 2*i;
//         bool current_test = (a_psum[i - 1] == expected);
//         if (!current_test)
//         {
//             printf("test %d failed. expected %d, actual %d\n", i, expected, a_psum[i]);
//         }
//         assert_filtered = assert_filtered && current_test;
//     }
//     printf("assert_filter=%s\n", assert_filtered ? "True" : "False");
//     free(a_psum);
//     free(a);
// }

int main()
{
    data_t *a = malloc(sizeof(data_t) * N);
    for (long  i = 0; i < N; i++)
        a[i] = i+1;
    data_t *a_psum = prefix_sum(a, N);
    filter_ret_t a_filtered = filter(a,N,even);
    data_t *a_filterd_array = a_filtered.filtered_array;
    index_t filtered_length = a_filtered.filtered_array_len;
    // for(index_t i=0;i<filtered_length;i++)
    //     printf("%d,",a_filterd_array[i]);
    // printf("\n");
    free(a_psum);
    free(a);
    free(a_filterd_array);
}