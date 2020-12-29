import math

class Function:
    def __init__(self, function):
        self.counter = 0
        self.f = function

    def __call__(self, x):
        self.counter += 1
        return self.f(x)


class Optimizer:
    def __init__(self):
        pass

    def optimize(self, function, a, b, derivative=None):
        x_opt = a
        f_opt = function(x_opt)
        """Здесь пишем код метода оптимизации"""        
        return x_opt, f_opt

class BruteOptimizer(Optimizer):
    def optimize(self, function, a, b, derivative=None):
        steps_num = 3000 # arbitrary

        x_opt = a
        f_opt = function(x_opt)

        x = a
        h = (b - a) / steps_num
        for i in range(0, steps_num):
            x += h
            y = function(x)
            if y < f_opt:
                x_opt = x
                f_opt = function(x_opt)

        return x_opt, f_opt


class GoldenSectionOptimizer(Optimizer):
    def optimize(self, function, a, b, derivative=None):
        delta = 1e-9 # arbitrary

        invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
        invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2

        h = b - a
        if h <= delta:
            x_opt = (a + b)/2
            f_opt = function(x_opt)
            return x_opt, f_opt

        # Number of steps required to reduce the search down to an delta-sized interval
        n = int(math.ceil(math.log(delta / h) / math.log(invphi)))

        c = a + invphi2 * h
        d = a + invphi * h
        yc = function(c)
        yd = function(d)

        for k in range(n-1):
            if yc < yd:
                b = d
                d = c
                yd = yc
                h = invphi * h
                c = a + invphi2 * h
                yc = function(c)
            else:
                a = c
                c = d
                yc = yd
                h = invphi * h
                d = a + invphi * h
                yd = function(d)

        if yc < yd:
            x_opt = (a + d)/2
        else:
            x_opt = (c + b)/2
        f_opt = function(x_opt)
        return x_opt, f_opt


class DichotomyOptimizer(Optimizer):
    def optimize(self, function, a, b, derivative=None):        
        delta = 1e-9

        mid = (a + b) / 2
        while (b - a > delta):

            if function(mid - delta) < function(mid + delta):
                b = mid
                mid = (a + b) / 2
            else:
                a = mid
                mid = (a + b) / 2

        x_opt = mid
        f_opt = function(x_opt)
        return x_opt, f_opt


class NewtonOptimizer(Optimizer):
    def optimize(self, function, a, b, derivative=None):
        x = (a + b)/2
        dfx = derivative(x)
        
        while a < x and x < b:
            dfx = derivative(x)
            if dfx == 0:
                break
            x = x - function(x)/dfx

        x_opt = x
        f_opt = function(x_opt)
        return x_opt, f_opt

def f1(x):
    return 2

def f2(x):
    return 10*x + 3

def f3(x):
    return 3*(x**2) - x + 2

def f4(x):
    return ((x - 3)**2 - 4)**2

def f5(x):
    return ((x - 3)**2 - 4)**2 + x

def f6(x):
    return math.sqrt(abs(x))


def main():
    def run_for_all_functions(optimizer_name, optimizer, a, b):
        print(f'--{optimizer_name}')
        run_and_output(optimizer, Function(f1), a, b)
        run_and_output(optimizer, Function(f2), a, b)
        run_and_output(optimizer, Function(f3), a, b)
        run_and_output(optimizer, Function(f4), a, b)
        run_and_output(optimizer, Function(f5), a, b)
        run_and_output(optimizer, Function(f6), a, b)

    def run_for_all_functions_with_derivatives(optimizer_name, optimizer, a, b):
        print(f'--{optimizer_name}')
        run_and_output(optimizer, Function(f1, df1), a, b)
        run_and_output(optimizer, Function(f2, df2), a, b)
        run_and_output(optimizer, Function(f3, df3), a, b)
        run_and_output(optimizer, Function(f4, df4), a, b)
        run_and_output(optimizer, Function(f5, df5), a, b)
        run_and_output(optimizer, Function(f6, df6), a, b)

    def run_and_output(optimizer, function, a, b):
        x, f = optimizer.optimize(function, a, b)
        print(f'The minimal value for {function.__name__} is {f} and is found at x = {x}')
        print(f'Function {function.__name__} was called {function.counter} times')
        print('--')

    brute_opt = BruteOptimizer()
    dichotomy_opt = DichotomyOptimizer()
    golden_section_opt = GoldenSectionOptimizer()
    newton_opt = NewtonOptimizer()

    a = 0
    b = 10

    run_for_all_functions('Brute force method', brute_opt, a, b)
    run_for_all_functions('Dichotomy method', dichotomy_opt, a, b)
    run_for_all_functions('Golden section method', golden_section_opt, a, b)
    run_for_all_functions('Newton method', newton_opt, a, b)
    

if __name__ == '__main__':
    main()
