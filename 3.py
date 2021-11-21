import random

# Генерация кривой
def GenEC():
    # характеристика поля
    p = random.randint(10**2, 10**3)
    while (true):
        p = random.randint(10**2, 10**3)
        if p in Primes():
            break
    # генерация конечного поля характеристики p
    K = GF(p)
    # выбор коэффициэнтов кривой
    A = K(randint(1,p))
    A = 4
    B = K(randint(1,p))
    B = 10
    # генерация кривой над конечным полем K
    E = EllipticCurve(K, [A,B])
    return E, p

def Schoof(E, p):
    # поиск простых чисел l_i
    ti = []
    li, m = ListLi(p)
    print('l_i: {}'.format(li))
    # вычисление значений для сравнений по модулю l_i
    for i in li:
        f = TraceOfFrobeniusModuloL(E, i, p)
        ti.append(f)
    print('T_i: {}'.format(ti))
    # китайская теорема об остатках 
    T = crt(ti, li)
    return T, m

# Поиск простых чисел Li
def ListLi(p):
    li = []
    m = 1
    fs = 4*sqrt(p)
    for i in Primes():
        m = m * i
        # из теоремы Хассе
        if m <= fs:
            li.append(i)
        else:
            break
    li.append(i)
    return li, m  

# след Фробениуса по mod l
def TraceOfFrobeniusModuloL(E, l, p):
    # построение расширения GF(p) = E.base_field()
    PR.<t> = PolynomialRing(E.base_field())
    # коэффициенты кривой
    A = E.a4()
    B = E.a6()
    
    if l != 2:
        # полином деления для l
        h = E.division_polynomial(l, t)
        # построение факторкольца F_p[x]/(h(x)) по идеалу полинома деления h(x)
        K.<x> = QuotientRing(PR, ideal(h)) 
        f = x^3 + A*x + B
        
        # (x^(p^2), y^(p^2)) + p(mod l) * (x, y) = T(mod l) * (x^p, y^p)
        # (x^p, y^p)
        xp = K(x**p)
        yp = K(f**((p-1) // 2))
        # (x^(p^2), y^(p^2))
        xp2, yp2 = xpypExp(xp, yp, K)     
        # вычисляем точку p_l
        pl_x, pl_y = Multiply(K(x), K(1), p, l, K, A, f)
        # вычисляем g = 𝜋^(2) + 𝑝 
        g_x, g_y = Summar(xp2, yp2, pl_x, pl_y, K, A, f)
        if g_y == 'Error':
            print('error h', l)
            h = gcd(h, lift(g_x))
            
            K.<x> = QuotientRing(PR, ideal(h))
            xp2, yp2 = xpypExp(xp, yp, K)
            # вычисляем точку p_l
            pl_x, pl_y = Multiply(K(x), K(1), p, l, K, A, f)
            # вычисляем g
            g_x, g_y = Summar(xp2, yp2, pl_x, pl_y, K, A, f)
            if g_y == 'Error':
                print('err h kk', l)

        # проверяем g
        if g_x == 0 or g_y == 0:
            return 0
        elif g_x == xp and g_y == yp:
            return 1
        elif g_x == (-1)*xp and g_y == (-1)*yp:
            return -1
        else:
            k = 2
            kxp, kyp = Multiply(xp, yp, k, l, K, A, f)
            if g_x == kxp and  g_y == kyp:
                    return k
              
            k += 1
            while k < l:
                kxp, kyp = Multiply(xp, yp, k, l, K, A, f)
                if g_x == kxp:
                    if g_y == kyp:
                        return k
                k += 1
        return 0
    else:
        f = t**3 + A*t + B
        # проверка является ли этот многочлен неприводимым над F2
        if f.is_irreducible():
            return 1
        else:
            return 0
    
# Вычисляем  квадрат эндоморфизма Фрабениуса 
def xpypExp(a, b, RR):
    xpp = Comp(a, a, RR)
    ypp = Comp(b, a, RR)
    ypp = RR(ypp * b)
    return xpp, ypp 

def Comp(a, b, RR):
    li = a.list()
    #print (li)
    res = 0
    d = 0
    for i in li:
        res = RR(res + i * (b^d))
        #print (res)
        d = d + 1
    return res

def Multiply(a, b, p, l, RR, A, f):
    d = p % l
    if d == 0:
        return 0, 0
    elif d == 1:
        return a, b
    elif d == 2:
        return Summar(a, b, a, b, RR, A, f)
    else:
        a0 = a
        b0 = b
        mask = bin(d)[3:]
        # print (mask)
        for i in range(len(mask)):
            a0, b0 = Summar(a0, b0, a0, b0, RR, A, f)
            if mask[i] == '1':
                a0, b0 = Summar(a0, b0, a, b, RR, A, f)
        return a0, b0
    
def Summar(a1, b1, a2, b2, RR, A, f):
    if a1 == a2 and b1 == b2:
        try:
            m = RR((3 * a1**2 + A) / (2 * b1 * f))
        except ZeroDivisionError:
            if RR(2 * b2 * f) == 0:
                return 1, 0  
            else:
                div = RR(2 * b1 * f)
                return div, 'Error'
    else:
        try:
            m = RR((b2 - b1) / (a2 - a1))
        except ZeroDivisionError:
            if RR(a2 - a1) == 0:
                return 1, 0
            div = RR(a2 - a1)
            return div, 'Error'

    a3 = RR(f * m**2 - a1 - a2)
    b3 = RR(m * (a1 - a3) - b1)
    return a3, b3   
    
def main():
    E, p = GenEC()
    print ('{}'.format(E))
    T, Z = Schoof(E, p)
    if T > (Z // 2):
        T = T - Z
    # число точек
    N = p + 1 - T 
    # print (E.trace_of_frobenius())
    # проверка равно ли полученное число точек порядку кривой
    if N == E.order():
        print ('Frobenius trace: T = {}'.format(T))
        print ('Order EC: N = {}'.format(N))
    else:
        print ('Error')

if __name__ == "__main__":
    main()