import subprocess
import tkinter as tk
from tkinter import messagebox

import matplotlib.pyplot as plt
import numpy as np

window = tk.Tk()

window_width = 1000
window_height = 800
x = (window.winfo_screenwidth() - window_width) // 2
y = (window.winfo_screenheight() - window_height) // 2

window.geometry(f"{window_width}x{window_height}+{x}+{y}")
window.resizable(False, False)
window.title("Задание параметров моделирования течения")
window.configure(bg='#fdfaff')


def process_data():
    messagebox.showinfo("Сообщение", "Программа выполняется. Необходимо подождать.")
    subprocess.run(['./alg.exe'] + [L.get(), R.get(), N1.get(), N2.get(), l1.get(), r.get(), l2.get(), l3.get(), Rk.get(), q.get(), Re.get(), eps.get(), delta.get(), it.get()])

    with open(r"data.txt", 'r') as file:
        file_L, file_R   = file.readline().split(' ')
        file_N1, file_N2 = file.readline().split(' ')
        file_l1, file_r, file_l2, file_l3, file_Rk = file.readline().split(' ')
        X = [float(number) for number in file.readline().split(' ')[:-1]]
        Y = [float(number) for number in file.readline().split(' ')[:-1]]

    with open(r"psi.txt", 'r') as file:
        psi = []
        for line in file:
            psi.append([float(number) for number in line.split(' ')[:-1]])

    file_L, file_R = float(file_L), float(file_R)
    file_N1, file_N2 = int(file_N1), int(file_N2)
    file_l1, file_r, file_l2, file_l3, file_Rk = float(file_l1), float(file_r), float(file_l2), float(file_l3), float(file_Rk) 
    x, y = np.array(X), np.array(Y)
    X, Y = np.meshgrid(x, y)
    psi = np.array(psi)

    plt.figure(figsize=(12, 5))
    plt.title('The flow in the blood vessel')
    plt.xlabel('length L')
    plt.ylabel('length R')

    # Верхнее сужение
    y1 = np.linspace(file_R - file_r, file_R + file_r, 100)
    x1, x2 = file_l1 - np.sqrt(file_r**2 - (y1 - file_R)**2), file_l1 + np.sqrt(file_r**2 - (y1 - file_R)**2)
    plt.fill_betweenx(y1, x1, x2, where=y1 <= file_R, facecolor='#090014', zorder=2)

    # Нижнее сужение
    y2 = np.linspace(-file_R - file_r, -file_R + file_r, 100)
    x1, x2 = file_l1 - np.sqrt(file_r**2 - (y2 + file_R)**2), file_l1 + np.sqrt(file_r**2 - (y2 + file_R)**2)
    plt.fill_betweenx(y2, x1, x2, where=y2 >= -file_R, facecolor='#090014', zorder=2)

    # Клапан
    if file_l2 and file_l3:
        plt.fill_betweenx(y, file_l2, file_l3, where=(y > -file_Rk) & (y < file_Rk), facecolor='#090014', zorder=2)

    # Линии функции тока
    plt.contour(X, Y, psi, 70, colors=['#3b0801', '#3b0801', '#4d0a02'], zorder=1)

    plt.show()


tk.Label(window, text="Параметры сосуда", width=20, anchor="w", padx=30, pady=15, font=("Arial", 20, "bold"), bg='#fdfaff').grid(row=0, column=0)

tk.Label(window, text="Длина L:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=1, column=0)
L = tk.Entry(window, width=7, font=("Arial", 14))
L.insert(0, '4,0')
L.grid(row=1, column=1, ipady=5)
tk.Label(window, text="Высота R:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=1, column=2)
R = tk.Entry(window, width=7, font=("Arial", 14))
R.insert(0, '1,0')
R.grid(row=1, column=3, ipady=5)

tk.Label(window, text="Число узлов N1:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=2, column=0)
N1 = tk.Entry(window, width=7, font=("Arial", 14))
N1.insert(0, 40)
N1.grid(row=2, column=1, ipady=5)
tk.Label(window, text="Число узлов N2:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=2, column=2)
N2 = tk.Entry(window, width=7, font=("Arial", 14))
N2.insert(0, 10)
N2.grid(row=2, column=3, ipady=5)

tk.Label(window, text="Параметры сужения", width=20, anchor="w", padx=30, pady=15, font=("Arial", 20, "bold"), bg='#fdfaff').grid(row=3, column=0)

tk.Label(window, text="Центр l1:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=4, column=0)
l1 = tk.Entry(window, width=7, font=("Arial", 14))
l1.insert(0, '1,0')
l1.grid(row=4, column=1, ipady=5)
tk.Label(window, text="Радиус:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=4, column=2)
r = tk.Entry(window, width=7, font=("Arial", 14))
r.insert(0, '0,4')
r.grid(row=4, column=3, ipady=5)

tk.Label(window, text="Параметры клапана", width=20, anchor="w", padx=30, pady=15, font=("Arial", 20, "bold"), bg='#fdfaff').grid(row=5, column=0)

tk.Label(window, text="Левая граница l2:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=6, column=0)
l2 = tk.Entry(window, width=7, font=("Arial", 14))
l2.insert(0, '1,7')
l2.grid(row=6, column=1, ipady=5)
tk.Label(window, text="Правая граница l3:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=7, column=0)
l3 = tk.Entry(window, width=7, font=("Arial", 14))
l3.insert(0, '2,0')
l3.grid(row=7, column=1, ipady=5)
tk.Label(window, text="Высота Rk:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=6, column=2)
Rk = tk.Entry(window, width=7, font=("Arial", 14))
Rk.insert(0, '0,7')
Rk.grid(row=6, column=3, ipady=5)

tk.Label(window, text="Параметры алгоритма", width=20, anchor="w", padx=30, pady=15, font=("Arial", 20, "bold"), bg='#fdfaff').grid(row=8, column=0)

tk.Label(window, text="Расход жидкости q:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=9, column=0)
q = tk.Entry(window, width=7, font=("Arial", 14))
q.insert(0, '1,0')
q.grid(row=9, column=1, ipady=5)
tk.Label(window, text="Число Рейнольдса Re:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=9, column=2)
Re = tk.Entry(window, width=7, font=("Arial", 14))
Re.insert(0, 20)
Re.grid(row=9, column=3, ipady=5)
tk.Label(window, text="Число eps:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=10, column=0)
eps = tk.Entry(window, width=7, font=("Arial", 14))
eps.insert(0, '0,00001')
eps.grid(row=10, column=1, ipady=5)
tk.Label(window, text="Число delta:", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=10, column=2)
delta = tk.Entry(window, width=7, font=("Arial", 14))
delta.insert(0, '0,0001')
delta.grid(row=10, column=3, ipady=5)
tk.Label(window, text="Число it (макс.):", width=20, anchor="center", padx=10, pady=10, font=("Arial", 16), bg='#fdfaff').grid(row=11, column=0)
it = tk.Entry(window, width=7, font=("Arial", 14))
it.insert(0, 3000)
it.grid(row=11, column=1, ipady=5)

submit_button = tk.Button(window, text="Построить график", width=20, pady=10, font=("Arial", 16), command=process_data)
submit_button.grid(row=12, column=0, columnspan=2, pady=40)

window.mainloop()
