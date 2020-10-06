# MPI
## Запуск программы
Перед компиляцией и запуском программы убедитесь, что установлена библиотека **MPI.h**. Установить MPI можно, например, с помощью команды:
```
$ sudo apt install mpich
```
### Компиляция
Для компиляции программы выполните команду:
```
$ make
```
После компиляции в директории появится исполняемый файл ***a.out***.
### Запуск
Запуск программы имеет вид:
```
$ mpirun -np p ./a.out n m
```
*n* - размер матрицы (целое положительное число)

*m* - размер блока (*m* должно быть меньше или равно *n*)

*p* - количество процессов (целое, большее 0)

## Результаты тестов
Контрольный запуск на 24 процессах, размере матрицы 25920 и размере блока 60 показал ускорение порядка ~21.5 раза.

----------------------------------------------
## Program run
Before compiling and running the program, make sure that the **MPI.h** library is installed. You can install MPI, for example, using the command:
### Compilation
To compile the program do:
```
$ make
```
After compilation, the ***a.out*** executable file will appear in the directory.
### Run
Program run looks like this:
```
$ mpirun -np p ./a.out n m
```
*n* - matrix size (positive integer)

*m* - block size (*m* should be greater or equal than *n*)

*p* - processes number (integer, greater than 0)

## Test results
A test run on 24 processes, 25920 matrix size and 60 block size showed an acceleration of the order of ~21.5 times.
