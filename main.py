import math


class SimplexMethod:

    # Сохраняет данные из файла, разбирая выражения на части
    def __init__(self, filename):
        file = open(filename, "r")
        self.mode = file.readline().strip()
        self.n = int(file.readline())
        self.f_coefs = file.readline().split(' ')
        self.m = int(file.readline())
        self.constraints_coefs = []
        self.constraints_signs = []
        self.constraints_limits = []
        for i in range(self.m):
            self.constraints_coefs.append([float(coef) for coef in file.readline().split(' ')])
            self.constraints_signs.append((int(file.readline())))
            self.constraints_limits.append(float(file.readline()))
            if self.constraints_limits[-1] < 0:
                self.constraints_limits[-1] *= -1
                for j in range(len(self.constraints_coefs[-1])):
                    self.constraints_coefs[-1][j] *= -1
                if self.constraints_signs[-1] == 0:
                    self.constraints_signs[-1] = 2
                else:
                    if self.constraints_signs[-1] == 2:
                        self.constraints_signs[-1] = 0
        self.f_coefs = [float(coef) for coef in self.f_coefs]
        self.table = []
        self.basic_vars = [0 for i in range(self.m)]

    def solve(self):
        self.form_start_table()
        self.do_iterations()
        self.print_result()

    # Формирует стартовый план
    def form_start_table(self):
        # Добавляем искусственные переменные
        for i in range(self.m):
            if self.constraints_signs[i] == 0:
                for j in range(self.m):
                    self.constraints_coefs[j].append(0)
                self.constraints_coefs[i][-1] = 1
            else:
                if self.constraints_signs[i] == 2:
                    for j in range(self.m):
                        self.constraints_coefs[j].append(0)
                    self.constraints_coefs[i][-1] = -1
                else:
                    for j in range(self.m):
                        self.constraints_coefs[j].append(0)
                    self.constraints_coefs[i][-1] = 1
        # Формируем строки таблицы, запоминая базисные переменные
        for i in range(self.m):
            for j in range(self.n, len(self.constraints_coefs[i])):
                if self.constraints_coefs[i][j] != 0:
                    row = [self.constraints_limits[i] / self.constraints_coefs[i][j]]
                    for l in range(len(self.constraints_coefs[i])):
                        row.append(self.constraints_coefs[i][l] / self.constraints_coefs[i][j])
                    self.table.append(row)
                    self.basic_vars[i] = j
        self.table.append([0])
        # Дописываем исходные коэффициенты функции
        if self.mode == "min":
            for coef in self.f_coefs:
                self.table[-1].append(-coef)
        else:
            for coef in self.f_coefs:
                self.table[-1].append(coef)
        while len(self.table[-1]) != len(self.table[0]):
            self.table[-1].append(0)
        # Дописываем коэффициенты функции с искусственными переменными
        check_row = [0 for j in range(len(self.table[0]))]
        if self.mode == "min":
            for i in range(len(self.table) - 1):
                for j in range(len(self.table[0])):
                    check_row[j] += self.table[i][j]
                    if j > self.n:
                        check_row[j] = 0
        else:
            for i in range(len(self.table) - 1):
                for j in range(len(self.table[0])):
                    check_row[j] -= self.table[i][j]
                    if j > self.n:
                        check_row[j] = 0
        self.table.append(check_row)
        print(f"Start table")
        for row in self.table:
            print(row)
        print("___________")
        return self.table

    # Итеративно оптимизирует план
    def do_iterations(self):
        iteration = 1
        while True:
            max = 0.0
            pivot_column = 0
            j = 1
            # Поиск ведущего столбца
            for j in range(1, self.n + 1):
                if self.table[-1][j] > max:
                    max = self.table[-1][j]
                    pivot_column = j
            # Проверка на оптимальность плана
            if abs(max) < 1e-10:
                for j in range(1, self.n + 1):
                    if self.table[-2][j] > max:
                        max = self.table[-2][j]
                        pivot_column = j
                if abs(max) < 1e-10:
                    return self.table
            # Поиск ведущей строки
            pivot_row = 0
            i = 0
            min_theta = math.inf
            for row in self.table:
                if i < self.m and row[pivot_column] > 0:
                    theta = row[0] / row[pivot_column]
                    if theta < min_theta:
                        min_theta = theta
                        pivot_row = i
                i += 1
            self.pivot(pivot_row, pivot_column)
            self.basic_vars[pivot_row] = pivot_column - 1
            print(f"Iteration: {iteration}")
            for row in self.table:
                print(row)
            print("___________")
            iteration += 1

    # Преобразование таблицы при переходе к следующему плану
    def pivot(self, row, column):
        j = 0
        # Все элементы ведущей строки делим на ведущий элемент
        pivot = self.table[row][column]
        for x in self.table[row]:
            self.table[row][j] = self.table[row][j] / pivot
            j += 1
        # В остальных строках:
        # элемент = элемент - коэффициент строки в ведущем столбце * новая ведущая строка в столбце элемента
        i = 0
        for xi in self.table:
            if i != row:
                ratio = xi[column]
                j = 0
                for xij in xi:
                    xij -= ratio * self.table[row][j]
                    self.table[i][j] = xij
                    j += 1
            i += 1

    # Вывод результата
    def print_result(self):
        if self.mode == "max":
            for j in range(len(self.table[-2])):
                self.table[-2][j] *= -1
        print("Final table:")
        for row in self.table:
            print(row)
        for i in range(len(self.basic_vars)):
            print(f"x{self.basic_vars[i] + 1} = {self.table[i][0]}")
        for i in range(self.n):
            if not self.basic_vars.__contains__(i):
                print(f"x{i + 1} = 0")
        print(f"Function extremum: {self.table[-2][0]}")


s = SimplexMethod('1 task.txt')
s.solve()
