# include <stdio.h>

int change(int *arr){
    arr[1] = 10;
    arr[4] = 8;
}

int main(){
    int arr[5] = {0};
    for (size_t i = 0; i < 5; i++)
    {
        printf("%d\n", arr[i]);
    }

    change(arr);

    for (size_t i = 0; i < 5; i++)
    {
        printf("%i\n", arr[i]);
    }

    return 0;
}