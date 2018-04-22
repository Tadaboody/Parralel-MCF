#include "omp.h"
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "pprefix.h"
typedef struct index_t_container
{
    index_t data;
} data_t;
typedef void* generic_p;
#define SIZE 10000000
index_t *prefix_sum(index_t *x, index_t n)
{
    //prefix sum happens _in place_. make sure not to free the array twice
    index_t *t = malloc(sizeof(data_t) * n);
    index_t i,j;
    for (j = 0; j < log2(n); j++) 
    {
    #pragma omp parallel private(i) //TODO: implement better
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

filter_ret_t filter(generic_p *array, index_t length, predicate p)
{
    index_t *bitmap, *bitsum;
    generic_p *filtered;
    bitmap = malloc(sizeof(index_t) * length);
    #pragma omp parralel for
    for(index_t i=0;i<length;i++)
    {
        bitmap[i] = p(array[i]);
    }
    bitsum = prefix_sum(bitmap, length);
    index_t filtered_length = bitsum[length - 1];
    filtered = malloc(sizeof(generic_p) * filtered_length);
    if(bitsum[0] > 0) // edge index
    {
        filtered[0] = array[0];
    }
    #pragma omp parallel for
    for (index_t i = 1; i < length; i++)
    {
        if(bitsum[i] > bitsum[i-1])
        {
            filtered[bitsum[i]-1] = array[i];
        }
    }
    free(bitsum);
    filter_ret_t ret = {filtered, filtered_length};
    return ret;
}

bool even(generic_p a)
{
    data_t* value = (data_t*) a;
    return (value->data) % 2 == 0;
}

void test_psum()
{
    bool assert_psum = true;
    index_t *a = malloc(sizeof(index_t) * SIZE);
    for (long i = 0; i < SIZE; i++)
        a[i] = i + 1;
    index_t *a_psum = prefix_sum(a, SIZE);
    #pragma omp parallel for reduction(&: assert_psum)
    for (index_t i = 1; i <= SIZE; i++)
    {
        index_t expected = (i * (i + 1)) / 2;
        bool current_test = (a_psum[i - 1] == expected);
        if (!current_test)
        {
            printf("test %d failed. expected %d, actual %d\n", i, expected, a_psum[i]);
        }
        assert_psum = assert_psum && current_test;
    }
    printf("assert_psum=%s\n", assert_psum ? "True" : "False");
    free(a_psum);
}

void test_filter()
{
    data_t **a = malloc(sizeof(data_t *) * SIZE);
    for (index_t i = 0; i < SIZE; i++)
    {
        a[i] = malloc(sizeof(data_t));
        a[i]->data = i + 1;
    }
    filter_ret_t a_filtered = filter((generic_p *)a, SIZE, even);
    data_t **a_filterd_array = a_filtered.filtered_array;
    index_t filtered_length = a_filtered.filtered_array_len;
    bool assert_filter = true;
    #pragma omp parallel for reduction(&& : assert_filter)
    for (index_t i = 0; i < filtered_length; i++)
    {
        data_t expected = {2 * (i+1)};
        data_t *actual = a_filterd_array[i];
        bool current_test = actual->data == expected.data;
        if (!current_test)
        {
            printf("test %d failed. expected %d, actual %d\n", i, expected.data, actual->data);
        }
        assert_filter = assert_filter && current_test;
    }
    printf("assert_filter=%s\n", assert_filter ? "True" : "False");
    free(a_filterd_array);
    for (index_t i = 0; i < SIZE; i++)
    {
        free(a[i]);
    }
    free(a);
}

// int main()
// {
//     test_psum();
//     test_filter();
// }