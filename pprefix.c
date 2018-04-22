#include "omp.h"
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
typedef long long index_t;
typedef struct index_t_container
{
    index_t data;
} data_t;
typedef void* generic_p;
#define N 10000000
typedef bool (*predicate)(generic_p);
typedef struct filter_ret_t
{
    generic_p filtered_array;
    index_t filtered_array_len;
}filter_ret_t;

index_t *prefix_sum(index_t *array, index_t n)
{
    index_t *x = malloc(sizeof(data_t) * n);
    index_t *t = malloc(sizeof(data_t) * n);
    index_t i,j;
    #pragma omp parallel for private(i)
    for (i = 0; i < n; i++)
        x[i] = array[i];
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

filter_ret_t filter(generic_p *array, generic_p *end,size_t jump, predicate p)
{
    index_t *bitmap, *bitsum;
    generic_p *filtered;
    size_t length = (end - array) / jump;
    size_t i;
    generic_p *member;
    bitmap = malloc(sizeof(index_t) * length);
    #pragma omp parralel for private(i,member)
    for (member = array, i = 0; member < end; member += jump, i++)
    {
        bitmap[i] = p(*member);
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
    free(bitmap);
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
    index_t *a = malloc(sizeof(index_t) * N);
    for (long i = 0; i < N; i++)
        a[i] = i + 1;
    index_t *a_psum = prefix_sum(a, N);
    #pragma omp parallel for reduction(&: assert_psum)
    for (index_t i = 1; i <= N; i++)
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
    free(a);
}

void test_filter()
{
    data_t **a = malloc(sizeof(data_t *) * N);
    for (index_t i = 0; i < N; i++)
    {
        a[i] = malloc(sizeof(data_t));
        a[i]->data = i + 1;
    }
    filter_ret_t a_filtered = filter((generic_p *)a, (generic_p *)(a + N), 1, even);
    data_t **a_filterd_array = a_filtered.filtered_array;
    index_t filtered_length = a_filtered.filtered_array_len;
    bool assert_filter = true;
    // #pragma omp parallel for reduction(&& : assert_filter)
    for (index_t i = 0; i < filtered_length; i++)
    {
        index_t expected = 2 * (i+1);
        index_t actual = a_filterd_array[i]->data;
        bool current_test = actual == expected;
        if (!current_test)
        {
            printf("test %d failed. expected %ld, actual %ld\n", i, expected, actual);
        }
        assert_filter = assert_filter && current_test;
    }
    printf("assert_filter=%s\n", assert_filter ? "True" : "False");
    free(a_filterd_array);
    for (index_t i = 0; i < N; i++)
    {
        free(a[i]);
    }
    free(a);
}
int main()
{
    test_psum();
    test_filter();
}