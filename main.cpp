#include <iostream>
#include <cstdlib>
#include <chrono>
using namespace std;
using namespace std::chrono;
void reset(int v[],int a[],int n)
{
    for(int i=0;i<n;i++)
    {
        v[i]=a[i];
    }
}
void printarray(int v[],int n)
{
    for(int i=0;i<n;i++)
    {
        cout<<v[i]<<" ";
    }
}
void bubbleSort(int v[],int n)
{
    int i,j;
    for (i=0;i<n-1;i++)
        for (j=0;j<n-i-1;j++)
            if (v[j]>v[j+1])
                swap(v[j],v[j+1]);
}
void insertionSort(int v[],int n)
{
    int i,key,j;
    for(i=1;i<n;i++)
    {
        key=v[i];
        j=i-1;
        while(j>=0 && v[j]>key)
        {
            v[j+1]=v[j];
            j=j-1;
        }
        v[j+1]=key;
    }
}
void selectionSort(int v[],int n)
{
    int i, j, min;
    for (i=0;i<n-1;i++)
    {
        min=i;
        for(j=i+1;j<n;j++)
        {
          if(v[j]<v[min])
              min=j;
        }
        if (min!=i)
            swap(v[min],v[i]);
    }
}
void merge(int array[], int const left, int const mid,int const right)
{
    auto const subArrayOne = mid - left + 1;
    auto const subArrayTwo = right - mid;
    auto *leftArray = new int[subArrayOne],
         *rightArray = new int[subArrayTwo];
    for (auto i = 0; i < subArrayOne; i++)
        leftArray[i] = array[left + i];
    for (auto j = 0; j < subArrayTwo; j++)
        rightArray[j] = array[mid + 1 + j];
    auto indexOfSubArrayOne
        = 0,
        indexOfSubArrayTwo
        = 0;
    int indexOfMergedArray
        = left;
    while (indexOfSubArrayOne < subArrayOne
           && indexOfSubArrayTwo < subArrayTwo) {
        if (leftArray[indexOfSubArrayOne]
            <= rightArray[indexOfSubArrayTwo]) {
            array[indexOfMergedArray]
                = leftArray[indexOfSubArrayOne];
            indexOfSubArrayOne++;
        }
        else {
            array[indexOfMergedArray]
                = rightArray[indexOfSubArrayTwo];
            indexOfSubArrayTwo++;
        }
        indexOfMergedArray++;
    }
    while (indexOfSubArrayOne < subArrayOne) {
        array[indexOfMergedArray]
            = leftArray[indexOfSubArrayOne];
        indexOfSubArrayOne++;
        indexOfMergedArray++;
    }
    while (indexOfSubArrayTwo < subArrayTwo) {
        array[indexOfMergedArray]
            = rightArray[indexOfSubArrayTwo];
        indexOfSubArrayTwo++;
        indexOfMergedArray++;
    }
    delete[] leftArray;
    delete[] rightArray;
}
void mergeSort(int array[], int const begin, int const end)
{
    if (begin >= end)
        return;
    auto mid = begin+(end - begin)/2;
    mergeSort(array,begin,mid);
    mergeSort(array,mid+1,end);
    merge(array,begin,mid,end);
}
int partition(int arr[], int low, int high)
{
    int pivot = arr[high];
    int i=(low-1);
    for (int j=low;j<=high-1;j++)
    {
        if (arr[j]<pivot)
        {
            i++;
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return (i + 1);
}
void quickSort(int arr[], int low, int high)
{
    if (low < high) {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}
void countSort(int A[],int SIZE)
{
  int HIGH=100001,C[HIGH];
  int B[SIZE + 1];
  for (int i = 0; i <= HIGH; i++)
    C[i] = 0;

  for(int j = 0; j < SIZE; j++)
    C[A[j]]++;

  for(int i = 1; i <= HIGH; i++)
    C[i] = C[i] + C[i - 1];

  for(int j = SIZE - 1; j >= 0; j--) {
    B[C[A[j]]] = A[j];
    C[A[j]]--;
  }

  for (int i = 1; i <= SIZE; i++)
    A[i - 1] = B[i];

}
int getMax(int arr[], int n)
{
    int mx = arr[0];
    for (int i = 1; i < n; i++)
        if (arr[i] > mx)
            mx = arr[i];
    return mx;
}
void countrSort(int arr[], int n, int exp)
{
    int output[n];
    int i, count[10] = { 0 };
    for (i = 0; i < n; i++)
        count[(arr[i] / exp) % 10]++;
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];
    for (i = n - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }
    for (i = 0; i < n; i++)
        arr[i] = output[i];
}
void radixSort(int arr[], int n)
{
    int m = getMax(arr, n);
    for (int exp = 1; m / exp > 0; exp *= 10)
        countrSort(arr, n, exp);
}
int shellSort(int arr[], int n)
{
    for (int gap = n/2; gap > 0; gap /= 2)
    {
        for (int i = gap; i < n; i += 1)
        {
            int temp = arr[i];
            int j;
            for (j = i; j >= gap && arr[j - gap] > temp; j -= gap)
                arr[j] = arr[j - gap];
            arr[j] = temp;
        }
    }
    return 0;
}
void heapify(int arr[], int N, int i)
{
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;
    if (l < N && arr[l] > arr[largest])
        largest = l;
    if (r < N && arr[r] > arr[largest])
        largest = r;
    if (largest != i) {
        swap(arr[i], arr[largest]);
        heapify(arr, N, largest);
    }
}
void heapSort(int arr[], int N)
{
    for (int i = N / 2 - 1; i >= 0; i--)
        heapify(arr, N, i);
    for (int i = N - 1; i > 0; i--)
    {
        swap(arr[0], arr[i]);
        heapify(arr, i, 0);
    }
}
int main()
{
    int n,i;
    cin>>n;
    int* v = (int*)malloc(n * sizeof(int)), *a = (int*)malloc(n * sizeof(int));
    for(i=0;i<n;i++)
        v[i]=rand() % 100000;
        reset(a,v,n);
    auto start = high_resolution_clock::now();
    bubbleSort(v,n);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout<<"bubble sort ="<<duration.count()<<endl;
    reset(v,a,n);
     start = high_resolution_clock::now();
    insertionSort(v,n);
     stop = high_resolution_clock::now();
     duration = duration_cast<microseconds>(stop - start);
    cout<<"insertion sort ="<<duration.count()<<endl;
    reset(v,a,n);
     start = high_resolution_clock::now();
    selectionSort(v,n);
     stop = high_resolution_clock::now();
     duration = duration_cast<microseconds>(stop - start);
    cout<<"selection sort ="<<duration.count()<<endl;
    reset(v,a,n);
     start = high_resolution_clock::now();
    mergeSort(v,0,n-1);
     stop = high_resolution_clock::now();
     duration = duration_cast<microseconds>(stop - start);
    cout<<"merge sort ="<<duration.count()<<endl;
    reset(v,a,n);
    start = high_resolution_clock::now();
    quickSort(v,0,n-1);
     stop = high_resolution_clock::now();
     duration = duration_cast<microseconds>(stop - start);
    cout<<"quick sort ="<<duration.count()<<endl;
    reset(v,a,n);
    start = high_resolution_clock::now();
    countSort(v,n);
     stop = high_resolution_clock::now();
     duration = duration_cast<microseconds>(stop - start);
    cout<<"count sort ="<<duration.count()<<endl;
    reset(v,a,n);
    start = high_resolution_clock::now();
    radixSort(v,n);
     stop = high_resolution_clock::now();
     duration = duration_cast<microseconds>(stop - start);
    cout<<"radix sort ="<<duration.count()<<endl;
    reset(v,a,n);
    start = high_resolution_clock::now();
    shellSort(v,n);
     stop = high_resolution_clock::now();
     duration = duration_cast<microseconds>(stop - start);
    cout<<"shell sort ="<<duration.count()<<endl;
    reset(v,a,n);
    start = high_resolution_clock::now();
    heapSort(v,n);
     stop = high_resolution_clock::now();
     duration = duration_cast<microseconds>(stop - start);
    cout<<"heap sort ="<<duration.count()<<endl;
    reset(v,a,n);
}
