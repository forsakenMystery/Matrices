import math as m


class Matrix:

    def __init__(self, lst):
        self.mat = lst
        self.n = len(lst)
        self.shape = [self.n]
        flag = True
        t = self.mat[0]
        while flag:
            try:
                s = len(t)
            except TypeError:
                break
            if s == 0:
                flag = False
            else:
                self.shape.append(s)
                t = t[0]

    def zeros(shape):
        a = [0]*shape[-1]
        for i in range(len(shape)-2, -1, -1):
            a = [a.copy() for i in range(shape[i])]
        return a

    def eye(shape):
        a = Matrix.zeros(shape)
        for i in range(len(a)):
            a[i][i] = 1
        return a

    def dot(x, y):
        a = len(x)
        b = len(y)
        c = len(y[0])
        z = Matrix.zeros([a, c])
        for i in range(a):
            for j in range(c):
                for k in range(b):
                    z[i][j] += x[i][k] * y[k][j]
        return z

    def LU_decomposition_pivoting(self):
        n, m = self.shape[0], self.shape[1]
        P = Matrix.eye([n, n])
        L = Matrix.eye([n, n])
        U = Matrix(self.mat).mat
        PF = Matrix.eye([n, n])
        LF = Matrix.zeros([n, n])
        for k in range(n - 1):
            mx = abs(U[k][k])
            index = k
            for i in range(k, n-1):
                if mx < abs(U[i][k]):
                    mx = abs(U[i][k])
                    index = i
            index = index + k
            if index != k:
                P = Matrix.eye([n, n])
                P[index][k:n], P[k][k:n] = P[k][k:n], P[index][k:n]
                U[index][k:n], U[k][k:n] = U[k][k:n], U[index][k:n]
                PF = Matrix.dot(P, PF)
                LF = Matrix.dot(P, LF)
            L = Matrix.eye([n, n])
            for j in range(k + 1, n):
                L[j][k] = -(U[j][k] / U[k][k])
                LF[j][k] = (U[j][k] / U[k][k])
            U = Matrix.dot(L, U)
        for i in range(len(LF)):
            LF[i][i] = 1
        return PF, LF, U

    def gauss(self, b):
        for y in range(self.shape[0]):
            maxrow = y
            for y2 in range(y+1, self.shape[0]):
                if abs(self.mat[y2][y])>abs(self.mat[maxrow][y]):
                    maxrow = y2
            self.mat[y], self.mat[maxrow] = self.mat[maxrow], self.mat[y]
            b[y], b[maxrow] = b[maxrow], b[y]
            print(self.mat)
            print(b)
            print("=====================")
            if abs(self.mat[y][y]) <= 1.e-6:
                return False, b
            for y2 in range(y+1, self.shape[0]):
                c = self.mat[y2][y] / self.mat[y][y]
                b[y2] -= b[y] * c
                for x in range(y, self.shape[1]):
                    self.mat[y2][x] -= self.mat[y][x] * c
        print("before")
        print(self.mat)
        print(b)
        print("=====================")
        for y in range(self.shape[0]-1, -1, -1):
            c = self.mat[y][y]
            print(self.mat)
            print(b)
            print("c", c)

            print("y", y)
            for y2 in range(y):
                print("()()()()()()()()()()()()()()()()")
                print("y2", y2)
                for x in range(self.shape[1]-1, y-1, -1):
                    print(self.mat)
                    print(b)
                    print("x", x)

                    print(self.mat[y2][y])
                    print(b[x] * self.mat[y2][y] / c)
                    b[y2] -= b[x] * self.mat[y2][x] / c
                    self.mat[y2][x] -= self.mat[y][x] * self.mat[y2][y] / c
                    print(self.mat)
                    print(b)
                    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7")
                print(self.mat)
                print(b)
            print("*****************")

            self.mat[y][y] /= c
            b[y] /= c

            for x in range(self.shape[0], self.shape[1]):
                self.mat[y][x] /= c

        return True, b

    def cholskey(self):
        result = Matrix.zeros(self.shape)
        print("result")
        print(result)
        for i in range(self.shape[0]):
            for j in range(i+1):
                res = 0
                for k in range(j):
                    res += result[i][k]*result[j][k]
                if i == j:
                    result[i][j] = m.sqrt(self.mat[i][i] - res)
                else:
                    result[i][j] = (1./result[j][j]*(self.mat[i][j] - res))
        print(result)

    def __repr__(self):
        return self.mat

    def __str__(self):
        return str(self.mat)


a = Matrix([[[2, 2], [3, 3]], [[2, 2], [3, 3]], [[2, 2], [3, 3]], [[2, 2], [3, 3]]])
print(a)
b = Matrix.zeros([2, 3])
print(b)
b[1][1] = 5
print(b)

print("cholskey")
A = Matrix([[6, 3, 4, 8], [3, 6, 5, 1], [4, 5, 10, 7], [8, 1, 7, 25]])
print(A.shape)
A.cholskey()

print("gauss")
A = Matrix([[0.5, 1.1, 3.1], [2.0, 4.5, 0.36], [5.0, 0.96, 6.5]])
A = Matrix([[1, 1, 1], [1, -1, 1], [2, 1, -1]])
solved, answer = A.gauss([4, 0, 12])
print(A)
if solved:
    print(answer)

A = [[1, 2], [3, 4]]
B = [[2, 0], [0, 2]]
print(Matrix.dot(A, B))
print("=FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF=")


A = [[2, 1, 1, 0], [4, 3, 3, 1], [8, 7, 9, 5], [6, 7, 9, 8]]
a = Matrix(A)
P1, _, _=a.LU_decomposition_pivoting()
print(P1)