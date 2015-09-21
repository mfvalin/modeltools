// An iterative implementation of quick sort
#include <stdio.h>
 
// A utility function to swap two elements
void swap ( int* a, int* b )
{
    int t = *a;
    *a = *b;
    *b = t;
//    printf("swap %d %d\n",*a,*b);
}

// A utility function to print contents of arr
void printArr( int arr[], int n )
{
    int i;
    for ( i = 0; i < n; ++i )
        printf( "%d ", arr[i] );
    printf("\n");
}
  
/* This function is same in both iterative and recursive*/
int partition (int arr[], int l, int h)
{
#if defined(ORIGINAL)
    int x = arr[h];
    int i = (l - 1);
    int j;
 
    for (j = l; j <= h- 1; j++)
    {
        if (arr[j] <= x)
        {
            i++;
            swap (&arr[i], &arr[j]);
        }
    }
    swap (&arr[i + 1], &arr[h]);
    return (i + 1);
#else
    int x = arr[(l+h)>>1];
//    int x = arr[l];
    int i = l;
    int j = h;
//    printf("partitioning %d to %d\n",l,h);
    while(1) {
      while(arr[j] > x) j--;
      while(arr[i] < x) i++;
//      printf("x=%d l=%d, h=%d, i=%d, j=%d\n",x,l,h,i,j);
      if(i<j) {
//        printf("S %d %d ",i,j);
        if(arr[i] == x && arr[j] == x) i++;
        else swap (&arr[i], &arr[j]);
      }
      else {
        printf("pivot = %d, pos = %d, len = %d\n",x,j,(h-l+1));
//        printArr(arr,33);
        return j;
      }
    }
#endif
}
 
/* A[] --> Array to be sorted, l  --> Starting index, h  --> Ending index */
void quickSortIterative (int arr[], int l, int h)
{
    // Create an auxiliary stack
    int stack[ h - l + 1 ];
 
    // initialize top of stack
    int top = -1;
    int maxdepth = 0;
    int pcalls = 0;
 
    // push initial values of l and h to stack
    stack[ ++top ] = l;
    stack[ ++top ] = h;
 
    // Keep popping from stack while is not empty
    while ( top >= 0 )
    {
        // Pop h and l
        h = stack[ top-- ];
        l = stack[ top-- ];
 
        // Set pivot element at its correct position in sorted array
        int p = partition( arr, l, h );
        pcalls++;
 
        // If there are elements on left side of pivot, then push left
        // side to stack
        if ( p-1 > l )
        {
            stack[ ++top ] = l;
            stack[ ++top ] = p - 1;
        }
 
        // If there are elements on right side of pivot, then push right
        // side to stack
        if ( p+1 < h )
        {
            stack[ ++top ] = p + 1;
            stack[ ++top ] = h;
        }
        maxdepth = top > maxdepth ? top : maxdepth;
    }
    printf("maxdepth = %d, partition calls = %d\n",(maxdepth+1)/2,pcalls);
}
 
// Driver program to test above functions
int main()
{
    int arr[] = {4, 3, 5, 2, 1, 3, 2, 3, 4, 3, 5, 2, 1, 3, 2, 3, 4, 3, 5, 2, 1, 3, 2, 3, 4, 3, 5, 2, 1, 3, 2, 3, 12};
//    int arr[64];
    int n = sizeof( arr ) / sizeof( *arr );
    int i;
//    for (i=0 ; i<n ; i++) arr[i] = n-i;
    printArr( arr, n ); printf("\n");
    quickSortIterative( arr, 0, n - 1 );
    printArr( arr, n ); printf("\n");
    return 0;
}
