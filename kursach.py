import numpy as np
%matplotlib inline
import matplotlib.pyplot as plt
%matplotlib inline
%config InlineBackend.figure_format = 'svg'
# Исходные данные
f1 = 60 # Частота 1-й гармоники 60
f2 = 80 # Частота 2-й гармоники 80
f3 = 100 # Частота 3-й гармоники 100
f4 = 120 # Частота 4-й гармоники 120
f5 = 140 # Частота 5-й гармоники 140
M = 5 # Кол-во гармоник 5
A = 2.5 # Амплитуда сигнала 2.5
N = 1024# Кол-во отсчетов
T = 0.001 # Период дискретизации 0.001
i = np.arange(N) # Временные отсчеты
X =(np.sin(2*np.pi*f1*T*i)+
    np.sin(2*np.pi*f2*T*i)+
    np.sin(2*np.pi*f3*T*i)+
    np.sin(2*np.pi*f4*T*i)+
    np.sin(2*np.pi*f5*T*i))*A/M
r = 12 # Разрядность 12
f = np.array([f1, f2, f3, f4, f5])

plt.figure(figsize=(10,6), dpi=300)
plt.plot(i*T, X)
plt.plot([i[0]*T, i[-1]*T], [A, A], "r--",alpha=0.5)
plt.plot([i[0]*T, i[-1]*T], [-A, -A], "r--",alpha=0.5)
plt.title('Полигармонический сигнал X(t)', fontsize=14)
plt.xlabel('t [сек]', fontsize=12)
plt.ylabel('X(t) [B]', fontsize=12)
plt.grid(True)
plt.show

# Генерируем шум
noise = np.random.normal(0, A*0.2, N)
# Строим график шума
plt.figure(figsize=(10,6), dpi=300)
plt.plot(i*T, noise)
plt.title('Сигнал шума noise(t)', fontsize=14)
plt.xlabel('t, [сек]', fontsize=12)
plt.ylabel('noise(t) [B]', fontsize=12)
plt.grid(True)
plt.show

# Накладываем шум на полезный сигнал и строим график
Xn = X + noise
fig, ax = plt.subplots(figsize=(10,6), dpi=300)
ax.set_title('Полигармонический сигнал X(t)', fontsize=14)
ax.set_xlabel('t [сек]', fontsize=12)
ax.set_xlabel('X(t) Xn(t) [B]', fontsize=12)
ax.plot(i*T, Xn, label='Сигнал с шумом Xn(t)')
ax.plot(i*T, X, label='Сигнал без шума X(t)')
ax.plot([i[0]*T,i[-1]*T], [A, A], 'r-*', alpha = 0.4)
ax.plot([i[0]*T,i[-1]*T], [-A, -A], 'r-*', alpha = 0.4)
ax.grid(True)
ax.legend()

out_file_name = 'output1.txt' # Имя входного файла
out_file = open(out_file_name, 'w') # Открытuе файла

# Вывод в файл
for p in range(N):
    out_file.write(f'{X[p]}\n')
# Поиск минимума и максимума
out_file.close()

# Модель аналого-цифрового преобразователя
Xmin = A # Верхняя граница X(t)
Xmax = -A # Нижняя граница X(t)
Xadc_max = np.power(2, r) - 1 # Верхняя граница АЦП
Xadc_min = 0 # Нижняя граница АЦП
# Чистый сигнал
C = Xadc_max-(X-Xmin)*(Xadc_max-Xadc_min)/(Xmax-Xmin)
c = C.round() # Округляем до целого числа
# Сигнал шумом
Cn = Xadc_max-(Xn-Xmin)*(Xadc_max-Xadc_min)/(Xmax-Xmin)
cn = Cn.round() # Округляем до целого числа
# Строим график оцифрованного сигнала
fig, ax = plt.subplots(figsize=(10,6),dpi=300)
ax.set_title('Сигнал АЦП ADC[tk]', fontsize=14)
ax.set_xlabel('tk [сек]', fontsize=12)
ax.set_ylabel('ADC[tk], ADCn[tk] [отсч]',fontsize=12)
ax.plot(i*T, cn, label='Сигнал с шумом Xn(t)')
ax.plot(i*T, c, label='Сигнал без шума X(t)')
ax.plot([i[0]*T, i[-1]*T], [Xadc_max, Xadc_max], 'r-*', alpha=0.3)
ax.plot([i[0]*T, i[-1]*T], [-Xadc_min, -Xadc_min], 'r-*', alpha=0.3)
ax.grid(True)
ax.legend()

# Обрезаем сигнал шума так, чтобы значения попадали в диапазон
cn[cn > Xadc_max] = Xadc_max
cn[cn < Xadc_min] = Xadc_min
fig, ax = plt.subplots(figsize=(10,6),dpi=300)
ax.set_title('Сигнал АЦП ADC[tk]', fontsize=14)
ax.set_xlabel('tk [сек]', fontsize=12)
ax.set_ylabel('ADC[tk], ADCn[tk] [отсч]',fontsize=12)
ax.plot(i*T, cn, label='Сигнал с шумом Xn(t)')
ax.plot(i*T, c, label='Сигнал без шума X(t)')
ax.plot([i[0]*T, i[-1]*T], [Xadc_max, Xadc_max], 'r-*', alpha=0.3)
ax.plot([i[0]*T, i[-1]*T], [-Xadc_min, -Xadc_min], 'r-*', alpha=0.3)
ax.grid(True)
ax.legend()

# Преобразование Фурье (приямое)
Xf=np.fft.fft(cn)/len(cn)
# Амплитудный спектр
Af=np.abs(Xf)
# Разрешение по частоте
b=1/N/T
# Массив частот
Ff = np.arange(N)*b

# График амплитудного спектра
plt.plot(Ff, Af)
plt.grid(True)
plt.show() 

# Поскольку на графике присутствует постоянная составляющая
# (Наблюдается всплеск на нулевой частоте)
# Рекомендуется выполнить центрирование ряда и повторить расчёт
cn_mean = np.mean(cn)
cn = cn-cn_mean

# Преобразование Фурье (прямое)
Xf=np.fft.fft(cn)/len(cn)
# Амлитудный спектр
Af=np.abs(Xf)

plt.plot(Ff, Af)
plt.grid(True)
plt.show()

plt.plot(Ff[:300], Af[:300])
plt.grid(True)
plt.show()

print(cn)
plt.plot(Ff[:300], Af[:300])
plt.grid(True)
plt.show()

# Определить уровень отсечки (фильтрации)
L = 50
annotation = 'L = '+str(L)

fig, ax = plt.subplots()
ax.plot(Ff, Af)
ax.plot([Ff[0], Ff[-1]], [L, L], 'r-', alpha=0.5)

ax.annotate(annotation, xy=(f[f.size-1]*10, L), xytext=(f[f.size-1]*10, L*1.5),
             ha='center', arrowprops={'arrowstyle':'->'})
plt.grid(True)
plt.show()

Wf = np.where(Af>L,1,0)
plt.plot(Ff, Wf)
plt.grid(True)

plt.plot(Ff[:150], Wf[:150])
plt.grid(True)

# Отрисовкаа фигуры в двух осях
plt.figure(figsize=(10, 6), dpi=300)
# Амлитудный спектр
plt.subplot(2,1,1)
plt.plot(Ff[:150], Af[:150])
plt.grid(True)
# Окно фильтра
plt.subplot(2,1,2)
plt.step(Ff[:150], Wf[:150])
plt.grid(True)

Np = 150
fig, ax = plt.subplots(figsize=(10,6), dpi = 300)
ax.plot(Ff[:Np], Af[:Np])
ax.plot([Ff[0], Ff[Np-1]], [L, L], 'r-', alpha = 0.5)
ax.step(Ff[:Np], Wf[:Np]*Af.max(), 'g-', alpha = 0.5)

ax.annotate(annotation, xy=(f[f.size-1]+f[0], L),
           xytext=(f[f.size-1]+f[0], L*1.5),
           ha='center', arrowprops={'arrowstyle':'->'})
plt.grid(True)
plt.show()

# Удаление составляющих комплесного спектра, соответствубщих шуму
Xf0 = Xf*Wf #Накладываем окно на спектр сигнала
# Выполняем обратное преобразование
# и оставляет только дейтсвительную часть
c0 = np.real(np.fft.ifft(Xf0))/T
# Восстанавливаем сигналы с учётои постоянной составляющей
c0 = c0 + cn_mean
cn = cn + cn_mean

# Визуализация
fig, ax = plt.subplots(figsize=(10,6), dpi=300)
ax.set_title('Зашумлённый и очищенный сигналы АЦП ADC[tk]')
ax.set_xlabel('tk [сек]', fontsize=12)
ax.set_ylabel('ADCn[tk], ADC0[tk] [отсч]', fontsize=12)
ax.plot(i*T, cn, label='Сигнал с шумом [tk]')
ax.plot(i*T, c0, label='Очищенный сигнал ADC0[tk]')
ax.plot([i[0]*T, i[-1]*T], [Xadc_max, Xadc_max], 'r-*', alpha=0.5)
ax.plot([i[0]*T, i[-1]*T], [-Xadc_min, -Xadc_min], 'r-*', alpha=0.5)
ax.grid(True)
ax.legend()

# Модель ЦАП
X0 = (Xadc_max - c0)*(Xmax-Xmin)/(Xadc_max-Xadc_min) + Xmin
# Визуализация
fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
ax.set_title('Исходный и востановленный сигнали X(t) X0(t)',fontsize=14)
ax.set_xlabel('t [сек]',fontsize=12)
ax.set_ylabel('X(t), X0(t) [B]', fontsize=12)
ax.plot(i*T, X, label='Исходный сигнал X(t)')
ax.plot(i*T, X0, label='Восстановленный сигнал X0(t)')
ax.plot([i[0]*T, i[-1]*T], [A, A], 'r-*', alpha=0.5)
ax.plot([i[0]*T, i[-1]*T], [-A, -A], 'r-*', alpha=0.5)
ax.grid(True)
ax.legend()

# Визуализация
fig, ax = plt.subplots(figsize=(10,6), dpi=300)
ax.set_title('Зашумленный и востановленный сигнали X(t) X0(t)',fontsize=14)
ax.set_xlabel('t [сек]',fontsize=12)
ax.set_ylabel('X(t), X0(t) [B]', fontsize=12)
ax.plot(i*T, Xn, label='Зашумленный сигнал Xn(t)')
ax.plot(i*T, X0, label='Восстановленный сигнал X0(t)')
ax.plot([i[0]*T, i[-1]*T], [A, A], 'r-*', alpha=0.5)
ax.plot([i[0]*T, i[-1]*T], [-A, -A], 'r-*', alpha=0.5)
ax.grid(True)
ax.legend()

# Расчёт и анализ ошибки
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
# Ошибка между сигнал с шумом и сигналом без шума
errorX = X - Xn
# Шибка между сигналом без шума и восстановленным сигналом
errorX0 = X - X0

# Визуализация
fig, ax = plt.subplots(figsize=(10,6), dpi=300)
ax.set_title('Ошибка между сигналами E(t)',fontsize=14)
ax.set_xlabel('t [сек]',fontsize=12)
ax.set_ylabel('E(t), X0(t) [B]', fontsize=12)
ax.plot(i*T, errorX, label='Ошибка между сигналом без шума и сигналом с шумом')
ax.plot(i*T, errorX0, label='Ошибка между сигналом без шума и аосстановленным')
ax.grid(True)
ax.legend()

# Расчёт статистических характеристик
print('Уровень L: {0}\n'.format(L))
print('MAE error X: {0}\n'.format(mean_absolute_error(X, Xn)))
print('MAE error X0: {0}\n'.format(mean_absolute_error(X, X0)))
print('MSE error X: {0}\n'.format(mean_squared_error(X, Xn)))
print('MSE error X0: {0}\n'.format(mean_squared_error(X, X0)))
print('R2 errorX: {0}\n'.format(r2_score(X, Xn)))
print('R2 errorX0: {0}\n'.format(r2_score(X, X0)))

# Визуализация
fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
ax.set_title('Исходный и востановленный сигналы X(t) X0(t)',fontsize=14)
ax.set_xlabel('t [сек]',fontsize=12)
ax.set_ylabel('X(t), X0(t) [B]', fontsize=12)
ax.plot(i[400:600]*T, X[400:600], label='Исходный сигнал X(t)')
ax.plot(i[400:600]*T, X0[400:600], label='Восстановленный сигнал X0(t)')
ax.grid(True)
ax.legend()

# Визуализация
fig, ax = plt.subplots(figsize=(10,6), dpi=300)
ax.set_title('Зашумленный и востановленный сигналы X(t) X0(t)',fontsize=14)
ax.set_xlabel('t [сек]',fontsize=12)
ax.set_ylabel('X(t), X0(t) [B]', fontsize=12)
ax.plot(i[400:600]*T, Xn[400:600], label='Зашумленный сигнал Xn(t)')
ax.plot(i[400:600]*T, X0[400:600], label='Восстановленный сигнал X0(t)')
ax.grid(True)
ax.legend()
