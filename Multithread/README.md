# Multithread
## Запуск программы
### Компиляция
Для компиляции программы выполните:
```
$ make
```
После компиляции в директории появится исполняемый файл ***a.out***.
### Запуск
Запуск программы имеет вид:
```
$ ./a.out n m p
```
*n* - размер матрицы (целое положительное число)

*m* - размер блока (*m* должно быть меньше или равно *n*)

*p* - количество потоков (целое, большее 0)
## Результаты тестов
Контрольный запуск на 24 потоках, размере матрицы 25920 и размере блока 60 показал ускорение порядка ~20.5 раз. 

------------------------------------------------
## Program run
### Compilation
To compile the program do:
```
$ make
```
After compilation, the ***a.out*** executable file will appear in the directory.
### Run
Program run looks like:
```
$ ./a.out n m p
```
*n* - matrix size (positive integer)

*m* - block size (*m* should be greater or equal than *n*)

*p* - threads number (integer, greater than 0)
## Test results
A test run on 24 threads, 25920 matrix size and 60 block size showed an acceleration of the order of ~20.5 times. 
