from time import process_time
import matplotlib.pyplot as plt
import numpy as np

import config as cnf


# второй множитель подинтегральной функции
# для вычисления образа - спектра сигнала
# - интеграл по delta t
def integ_direct(x_j, h, k):
    return 2.0 * h * (1 - np.cos(k*h)) / (k * h)**2 * np.exp(-1j * k * x_j)


# ...
# для вычисления прообраза - сигнала
# ...
def integ_inverse(x_j, h, k):
    return 2.0 * h * (1 - np.cos(k*h)) / (k * h)**2 * np.exp(1j * k * x_j)


# ===========================================================================
start = process_time()

t = []  # массив точек по оси времени t
p = []  # массив соответсвующих значений сигнала

with open(cnf.file, 'r') as data:
    counter = 0
    for line in data:
        if line[0] == '%':
            counter += 1
            continue
        items = line.split('\t')
        t.append(float(items[0]))
        p.append(float(items[1]))
        counter += 1

t = [i*(10**-3) for i in t]  # переводим микроСекунды (мкс) в миллиСекунды (мс)

# --- визуализация -------------------------------
fig, ax = plt.subplots(figsize=[12, 5])
ax.plot(t, p, label='исходный сигнал')
ax.set_xlabel('Время [мс]')
ax.set_ylabel('Напряжение [мВ]')
ax.set_title("График сигнала")
ax.grid(linestyle='-.')
ax.legend()
# plt.show()

# ----------------------- Спектр сигнала -------------------------------------
# прореживаем сигнал
t = [item for i, item in enumerate(t) if (i+1) % cnf.k_thin == 0]
p = [item for i, item in enumerate(p) if (i+1) % cnf.k_thin == 0]
th = t[1] - t[0]

# Вычисляем спектр через ППФ
f = np.arange(cnf.f_start, cnf.f_end, cnf.f_step)  # размерная частота [кГц]

omega = [i*2*np.pi for i in f]  # безразмерная частота
mgh = omega[1] - omega[0]  # шаг по omega

p_ft = []  # спектр сигнала, образ
for mg in omega:
    if cnf.mode == 'direct':
        p_ft.append(np.sum([p_j * np.exp(complex(0, -1) * mg * t_j)
                            for p_j, t_j in zip(p, t)]) * th)
    elif cnf.mode == 'splines':
        p_ft.append(np.sum([p_j * integ_direct(t_j, th, mg)
                            for p_j, t_j in zip(p, t)]))
p_ft_abs = list(map(abs, p_ft))

# --- визуализация -------------------------------
fig1, ax1 = plt.subplots(figsize=[12, 5])
markerline, stemline, baseline, = ax1.stem(f, p_ft_abs)
plt.setp(stemline, linewidth=1.25)
plt.setp(markerline, markersize=3)
ax1.set_xlabel('Частота [кГц]')
ax1.set_ylabel('Напряжение [мВ]')
ax1.set_title("График амплитудного спектра сигнала")
ax1.grid(linestyle='-.')
# plt.show()

# ----------------------- Сигнал после ОПФ -------------------------------------
p_ift = []
for t_j in t:
    if cnf.mode == 'direct':
        p_ift.append(np.sum([pw_i * np.exp(complex(0, 1) * mg_i * t_j)
                             for pw_i, mg_i in zip(p_ft, omega)]) * mgh)
    elif cnf.mode == 'splines':
        p_ift.append(np.sum([pw_i * integ_inverse(mg_i, mgh, t_j)
                             for pw_i, mg_i in zip(p_ft, omega)]))
p_ift_re = [i.real / np.pi for i in p_ift]

# --- визуализация -------------------------------
ax.plot(t, p_ift_re, label='сигнал после ОПФ')
ax.set_xlabel('Время [мс]')
ax.set_ylabel('Напряжение [мВ]')
ax.set_title("График сигнала")
ax.grid(linestyle='-.')
ax.legend()

end = process_time()
print(f'time: {end-start: .5f} s')

plt.show()
