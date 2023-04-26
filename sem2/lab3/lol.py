import numpy as np

elem = np.array([1, 2])
arr = np.array([elem])  # массив массивов

new_element = np.array([7, 8])  # новый элемент, который мы хотим добавить

arr = np.append(arr, [new_element])  # добавляем новый элемент в конец массива

print(arr)