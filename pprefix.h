#include <stdbool.h>
typedef long index_t;
typedef void* generic_p;
typedef bool (*predicate)(generic_p);
typedef struct filter_ret_t
{
    generic_p filtered_array;
    index_t filtered_array_len;
}filter_ret_t;
index_t *prefix_sum(index_t *x, index_t n);
filter_ret_t filter(generic_p *array, index_t length, predicate p);